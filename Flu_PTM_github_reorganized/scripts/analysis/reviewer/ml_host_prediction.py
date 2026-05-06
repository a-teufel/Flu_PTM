# Machine-learning host prediction from PTM states
# Teufel et al. (2025) Genome Biology and Evolution — reviewer revision
#
# For each protein, trains Random Forest and Logistic Regression classifiers
# to predict broad host category (Human / Avian / Swine) from the PTM symbol
# profile at all sites, fast-evolving sites only, and slow-evolving sites only.
# Uses ordinal-encoded PTM symbols as features; 5-fold stratified CV.
#
# Outputs written to results/reviewer_analyses/:
#   ml_results_summary.csv      — balanced accuracy per protein/model/feature set
#   ml_feature_importance.csv   — Random Forest site importances
#   ml_confusion_matrices.csv   — confusion matrix entries (long format)
#
# Usage (from repo root):
#   python scripts/analysis/ml_host_prediction.py

import re, json, csv, warnings
import numpy as np
import pandas as pd
from pathlib import Path
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import LabelEncoder, OrdinalEncoder, StandardScaler, OneHotEncoder
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.metrics import balanced_accuracy_score, accuracy_score, confusion_matrix
from sklearn.pipeline import Pipeline

warnings.filterwarnings("ignore")

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE   = Path(__file__).resolve().parents[2]
MULTI  = BASE / "data/ptm/multistate"
MAPS   = BASE / "data/ptm/mappings"
META_D = BASE / "data/raw_sequences"
PAPER_RATE_FILE = MAPS / "paper_site_rate_annotations.csv"
OUT    = BASE / "results/reviewer_analyses"
OUT.mkdir(parents=True, exist_ok=True)

PROTEINS = [
    "PB2_polymerase", "PB1_polymerase", "PB1-F2_protein",
    "PA_polymerase",  "PA-X_protein",   "NP_protein",
    "NA_neuraminidase","M1_matrix_protein","M2_ion_channel",
    "NS1_protein",    "NS2_protein",    "NEP_protein",    "HA_full",
]

AVIAN_RE = re.compile(
    r"Gallus|Anas|Anatidae|Aves|Meleagris|Numididae|Spatula|Mareca|Cygnus|"
    r"Columbidae|Phasian|Accipitriformes|Psittacidae|Arenaria|Cairina", re.I)

PTM_SYMBOLS = sorted("NQRSTXY-?")


def broad_host(h: str):
    if not h or h.strip() in ("", "Unknown"):
        return None
    if re.search(r"Homo sapiens|human", h, re.I):
        return "Human"
    if re.search(r"Sus scrofa|pig|swine|porcine", h, re.I):
        return "Swine"
    if AVIAN_RE.search(h):
        return "Avian"
    return None


def load_metadata():
    meta = {}
    for strain in ("H1N1", "H5N1", "H7N9"):
        path = META_D / strain / f"{strain}_seq_info.csv"
        with open(path) as f:
            for row in csv.DictReader(f):
                acc = row["Accession"].strip().rstrip(".1")
                meta[acc] = {"strain": strain,
                             "host_broad": broad_host(row.get("Host", ""))}
    return meta


def parse_multistate(path: Path):
    lines = open(path).readlines()
    start = next(i for i, l in enumerate(lines) if l.strip().startswith("MATRIX")) + 1
    end   = next(i for i, l in enumerate(lines[start:], start) if l.strip() == ";")
    out = {}
    for ln in lines[start:end]:
        parts = ln.split()
        if len(parts) >= 2:
            out[parts[0]] = parts[1:]
    return out


def get_char_labels(path: Path):
    lines = open(path).readlines()
    labels, in_block = [], False
    for ln in lines:
        t = ln.strip()
        if t.startswith("CHARSTATELABELS"):
            in_block = True
            continue
        if in_block:
            if t.startswith(";"):
                break
            m = re.search(r"(\d+)\s+pos_(\d+)", t)
            if m:
                labels.append(int(m.group(2)))
    return labels


def make_resolver(json_map: dict):
    def resolve(taxon: str) -> str:
        if taxon in json_map:
            v = json_map[taxon]
            m = re.match(r"^(H1N1|H5N1|H7N9)__([A-Z0-9]+)_$", v)
            return m.group(2) if m else re.sub(r"\.1$", "", v)
        m = re.match(r"^(H1N1|H5N1|H7N9)__([A-Z0-9]+)_$", taxon)
        return m.group(2) if m else re.sub(r"\.1$", "", taxon)
    return resolve


def load_paper_rate_annotations(path: Path):
    if not path.exists():
        raise FileNotFoundError(f"Missing canonical paper site file: {path}")

    ann = pd.read_csv(path)
    required = {"protein", "site_idx", "rate_class"}
    if not required.issubset(ann.columns):
        raise ValueError(
            f"{path} is missing required columns: {sorted(required)}"
        )

    site_map = {}
    for protein in PROTEINS:
        sub = ann[ann["protein"] == protein]
        fast = set(sub.loc[sub["rate_class"].str.lower() == "fast", "site_idx"].astype(int))
        slow = set(sub.loc[sub["rate_class"].str.lower() == "slow", "site_idx"].astype(int))
        site_map[protein] = {"fast": fast, "slow": slow}
    return site_map


PAPER_SITE_MAP = load_paper_rate_annotations(PAPER_RATE_FILE)


def load_fast_slow_sites(protein: str):
    ann = PAPER_SITE_MAP.get(protein, {"fast": set(), "slow": set()})
    return ann["fast"], ann["slow"]


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
meta = load_metadata()

results_rows    = []
importance_rows = []
cm_rows         = []

cv5 = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

for protein in PROTEINS:
    mfile = MULTI / f"{protein}_multistate.txt"
    if not mfile.exists():
        continue

    json_path = MAPS / f"{protein}_id_mapping.json"
    json_map  = json.load(open(json_path)) if json_path.exists() else {}
    resolve   = make_resolver(json_map)

    records     = parse_multistate(mfile)
    char_labels = get_char_labels(mfile)
    n_sites     = max(len(v) for v in records.values()) if records else 0
    fast_idx, slow_idx = load_fast_slow_sites(protein)

    # Build (X, y) for this protein
    rows_X, rows_y, rows_strain = [], [], []
    for taxon, syms in records.items():
        acc  = resolve(taxon)
        info = meta.get(acc, {})
        hb   = info.get("host_broad")
        st   = info.get("strain")
        if hb not in ("Human", "Avian", "Swine") or st not in ("H1N1", "H5N1", "H7N9"):
            continue
        feats = [syms[si] if si < len(syms) else "?" for si in range(n_sites)]
        rows_X.append(feats)
        rows_y.append(hb)
        rows_strain.append(st)

    if len(rows_X) < 30:
        print(f"  skip {protein} — only {len(rows_X)} labeled samples")
        continue

    X_raw = np.array(rows_X)
    y     = np.array(rows_y)
    strain_arr = np.array(rows_strain).reshape(-1, 1)

    # Require ≥ 10 samples per class for reliable CV
    vc = dict(zip(*np.unique(y, return_counts=True)))
    if min(vc.values()) < 10:
        print(f"  skip {protein} — class too small: {vc}")
        continue

    # Ordinal-encode PTM symbols (unknown → -1)
    oe    = OrdinalEncoder(handle_unknown="use_encoded_value",
                           unknown_value=-1,
                           categories=[PTM_SYMBOLS] * n_sites)
    X_enc = oe.fit_transform(X_raw)

    try:
        ohe = OneHotEncoder(handle_unknown="ignore", sparse_output=False)
    except TypeError:
        ohe = OneHotEncoder(handle_unknown="ignore", sparse=False)
    X_strain = ohe.fit_transform(strain_arr)

    le_host = LabelEncoder()
    y_enc   = le_host.fit_transform(y)
    classes = le_host.classes_

    fast_cols = [i for i in range(n_sites) if (i + 1) in fast_idx]
    slow_cols = [i for i in range(n_sites) if (i + 1) in slow_idx]
    feature_sets = {"all": list(range(n_sites))}
    if fast_cols:
        feature_sets["fast"] = fast_cols
    if slow_cols:
        feature_sets["slow"] = slow_cols

    models = {
        "RandomForest": RandomForestClassifier(
            n_estimators=200, max_depth=8,
            class_weight="balanced", random_state=666, n_jobs=1),
        "LogisticReg":  Pipeline([
            ("sc",  StandardScaler()),
            ("clf", LogisticRegression(
                C=0.5, max_iter=500, class_weight="balanced",
                solver="lbfgs", random_state=666))]),
    }

    # Baseline model using strain only
    for model_name, clf in models.items():
        try:
            y_pred = cross_val_predict(clf, X_strain, y_enc, cv=cv5)
            acc    = accuracy_score(y_enc, y_pred)
            bacc   = balanced_accuracy_score(y_enc, y_pred)
            cm     = confusion_matrix(y_enc, y_pred, labels=range(len(classes)))
            results_rows.append({
                "protein":           protein,
                "feature_set":       "strain_only",
                "model":             model_name,
                "n_samples":         X_strain.shape[0],
                "n_sites_used":      0,
                "n_strain_features": X_strain.shape[1],
                "accuracy":          round(acc, 4),
                "balanced_accuracy": round(bacc, 4),
            })
            for i, tc in enumerate(classes):
                for j, pc in enumerate(classes):
                    cm_rows.append({"protein": protein, "feature_set": "strain_only",
                                    "model": model_name, "true": tc,
                                    "predicted": pc, "count": int(cm[i, j])})
        except Exception as e:
            print(f"  WARN {protein}/strain_only/{model_name}: {e}")

    for feat_name, cols in feature_sets.items():
        X_ptm = X_enc[:, cols]
        X_adj = np.hstack([X_ptm, X_strain])
        for model_name, clf in models.items():
            # PTM-only model
            try:
                y_pred = cross_val_predict(clf, X_ptm, y_enc, cv=cv5)
                acc    = accuracy_score(y_enc, y_pred)
                bacc   = balanced_accuracy_score(y_enc, y_pred)
                cm     = confusion_matrix(y_enc, y_pred, labels=range(len(classes)))
                results_rows.append({
                    "protein":           protein,
                    "feature_set":       f"ptm_{feat_name}",
                    "model":             model_name,
                    "n_samples":         X_ptm.shape[0],
                    "n_sites_used":      X_ptm.shape[1],
                    "n_strain_features": 0,
                    "accuracy":          round(acc, 4),
                    "balanced_accuracy": round(bacc, 4),
                })
                for i, tc in enumerate(classes):
                    for j, pc in enumerate(classes):
                        cm_rows.append({"protein": protein, "feature_set": f"ptm_{feat_name}",
                                        "model": model_name, "true": tc,
                                        "predicted": pc, "count": int(cm[i, j])})
            except Exception as e:
                print(f"  WARN {protein}/ptm_{feat_name}/{model_name}: {e}")

            # Strain-adjusted model (PTM + strain covariates)
            try:
                y_pred = cross_val_predict(clf, X_adj, y_enc, cv=cv5)
                acc    = accuracy_score(y_enc, y_pred)
                bacc   = balanced_accuracy_score(y_enc, y_pred)
                cm     = confusion_matrix(y_enc, y_pred, labels=range(len(classes)))
                results_rows.append({
                    "protein":           protein,
                    "feature_set":       f"ptm_plus_strain_{feat_name}",
                    "model":             model_name,
                    "n_samples":         X_adj.shape[0],
                    "n_sites_used":      X_ptm.shape[1],
                    "n_strain_features": X_strain.shape[1],
                    "accuracy":          round(acc, 4),
                    "balanced_accuracy": round(bacc, 4),
                })
                for i, tc in enumerate(classes):
                    for j, pc in enumerate(classes):
                        cm_rows.append({"protein": protein, "feature_set": f"ptm_plus_strain_{feat_name}",
                                        "model": model_name, "true": tc,
                                        "predicted": pc, "count": int(cm[i, j])})
            except Exception as e:
                print(f"  WARN {protein}/ptm_plus_strain_{feat_name}/{model_name}: {e}")

    # RF feature importance — all PTM sites adjusted for strain covariates
    rf = RandomForestClassifier(n_estimators=300, max_depth=8,
                                 class_weight="balanced", random_state=42, n_jobs=1)
    X_rf = np.hstack([X_enc, X_strain])
    rf.fit(X_rf, y_enc)
    importances_ptm = rf.feature_importances_[:n_sites]
    for ci, imp in enumerate(importances_ptm):
        orig_pos = char_labels[ci] if ci < len(char_labels) else ci + 1
        importance_rows.append({
            "protein":  protein,
            "site_idx": ci + 1,
            "orig_pos": orig_pos,
            "importance": round(float(imp), 6),
            "importance_adjusted_for_strain": True,
            "is_fast":  (ci + 1) in fast_idx,
            "is_slow":  (ci + 1) in slow_idx,
        })

    vc_str = "  ".join(f"{k}:{v}" for k, v in vc.items())
    print(f"  {protein:25s}  {n_sites} sites  {vc_str}")

# ---------------------------------------------------------------------------
# Write outputs
# ---------------------------------------------------------------------------
pd.DataFrame(results_rows).to_csv(OUT / "ml_results_summary.csv", index=False)
pd.DataFrame(importance_rows).to_csv(OUT / "ml_feature_importance.csv", index=False)
pd.DataFrame(cm_rows).to_csv(OUT / "ml_confusion_matrices.csv", index=False)
print(f"\nSaved outputs to {OUT}")
