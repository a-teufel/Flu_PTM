# Machine-learning host prediction from PTM states
# Teufel et al. (2025) Genome Biology and Evolution — reviewer revision
#
# For each protein, trains Random Forest and Logistic Regression classifiers
# to predict broad host category (Human / Avian / Swine / Mammalian_Other)
# from the PTM symbol
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

# By default, run models on all PTM sites only.
# Set to True if you also want fast/slow subset models when annotations are available.
USE_FAST_SLOW_SUBSETS = True

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE   = Path(__file__).resolve().parents[3]
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

# Exact host-to-class mapping derived from every unique Host value in
# data/raw_sequences/{H1N1,H5N1,H7N9}/*_seq_info.csv.  Nothing is guessed.
HOST_MAP = {
    "Homo sapiens":           "Human",
    "Sus scrofa":             "Swine",
    # Avian
    "Gallus gallus":          "Avian",
    "Anatidae":               "Avian",
    "Aves":                   "Avian",
    "Meleagris gallopavo":    "Avian",
    "Anas platyrhynchos":     "Avian",
    "Phasianinae":            "Avian",
    "Mareca americana":       "Avian",
    "Numididae":              "Avian",
    "Anas carolinensis":      "Avian",
    "Spatula clypeata":       "Avian",
    "Cygnus olor":            "Avian",
    "Columbidae":             "Avian",
    "Arenaria interpres":     "Avian",
    "Cairina moschata":       "Avian",
    "Spatula discors":        "Avian",
    "Numididae sp.":          "Avian",
    "Anas cyanoptera":        "Avian",
    "Parus major":            "Avian",
    "Psittacidae":            "Avian",
    "Accipitriformes":        "Avian",
    "Passer montanus":        "Avian",
    "Accipitridae":           "Avian",
    "Tadorna":                "Avian",
    "Hirundo rustica":        "Avian",
    "Aythya ferina":          "Avian",
    "Galliformes":            "Avian",
    "Tachybaptus ruficollis": "Avian",
    # Mammalian_Other (non-human, non-swine mammals)
    "Panthera tigris":        "Mammalian_Other",
    "Mus musculus":           "Mammalian_Other",
    "Mustela lutreola":       "Mammalian_Other",
    "Mustela putorius furo":  "Mammalian_Other",
    "Suricata suricatta":     "Mammalian_Other",
    "Felis catus":            "Mammalian_Other",
    # Empty / unknown
    "":                       None,
}

HOST_CLASSES = ("Human", "Avian", "Swine", "Mammalian_Other")


def broad_host(h: str):
    """Map a raw Host string to a broad class using the exact lookup table.
    Returns None for empty/unknown hosts (they will be excluded from analysis)."""
    return HOST_MAP.get(h.strip() if h else "", None)


def load_metadata():
    meta = {}
    host_counts = {}
    host_to_class = {}
    for strain in ("H1N1", "H5N1", "H7N9"):
        path = META_D / strain / f"{strain}_seq_info.csv"
        with open(path) as f:
            for row in csv.DictReader(f):
                acc = row["Accession"].strip().rstrip(".1")
                raw_host = (row.get("Host", "") or "").strip()
                hb = broad_host(raw_host)
                meta[acc] = {"strain": strain, "host_broad": hb, "host_raw": raw_host}
                host_counts[raw_host] = host_counts.get(raw_host, 0) + 1
                if raw_host not in host_to_class:
                    host_to_class[raw_host] = hb

    host_map_rows = []
    for host, n in sorted(host_counts.items(), key=lambda x: (-x[1], x[0])):
        host_map_rows.append({
            "host_raw": host,
            "host_broad": host_to_class.get(host),
            "n_records": n,
        })
    return meta, host_map_rows


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


def read_ptm_symbols(path: Path):
    """Read the canonical PTM symbol alphabet from the SYMBOLS= line in the
    nexus file header.  This is the authoritative source — no inference needed."""
    with open(path) as f:
        for line in f:
            m = re.search(r'SYMBOLS="([^"]+)"', line)
            if m:
                symbols = sorted(m.group(1))  # individual characters
                if "?" not in symbols:
                    symbols.append("?")
                    symbols.sort()
                return symbols
    raise ValueError(f"No SYMBOLS= line found in {path}")


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
        print(f"WARN missing site-rate annotation file: {path} (fast/slow subsets disabled)")
        return {}

    ann = pd.read_csv(path)
    required = {"protein", "site_idx", "rate_class"}
    if not required.issubset(ann.columns):
        print(f"WARN annotation file missing required columns {sorted(required)}; fast/slow subsets disabled")
        return {}

    site_map = {}
    for protein in PROTEINS:
        sub = ann[ann["protein"] == protein]
        fast = set(sub.loc[sub["rate_class"].str.lower() == "fast", "site_idx"].astype(int))
        slow = set(sub.loc[sub["rate_class"].str.lower() == "slow", "site_idx"].astype(int))
        site_map[protein] = {"fast": fast, "slow": slow}
    return site_map


PAPER_SITE_MAP = load_paper_rate_annotations(PAPER_RATE_FILE)


def load_fast_slow_sites(protein: str):
    if not USE_FAST_SLOW_SUBSETS:
        return set(), set()
    ann = PAPER_SITE_MAP.get(protein, {"fast": set(), "slow": set()})
    return ann["fast"], ann["slow"]


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
meta, host_map_rows = load_metadata()

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
        if hb not in HOST_CLASSES or st not in ("H1N1", "H5N1", "H7N9"):
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

    # Keep only classes with enough support for stable CV within each protein.
    vc = dict(zip(*np.unique(y, return_counts=True)))
    eligible_classes = [k for k, v in vc.items() if v >= 10]
    if len(eligible_classes) < 2:
        print(f"  skip {protein} — insufficient class support: {vc}")
        continue

    keep = np.isin(y, eligible_classes)
    X_raw = X_raw[keep]
    y = y[keep]
    strain_arr = strain_arr[keep]
    vc = dict(zip(*np.unique(y, return_counts=True)))

    # Ordinal-encode PTM symbols using the canonical alphabet from the nexus header.
    ptm_symbols = read_ptm_symbols(mfile)
    oe    = OrdinalEncoder(handle_unknown="use_encoded_value",
                           unknown_value=-1,
                           categories=[ptm_symbols] * n_sites)
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
pd.DataFrame(host_map_rows).to_csv(OUT / "ml_host_mapping_audit.csv", index=False)
print(f"\nSaved outputs to {OUT}")
