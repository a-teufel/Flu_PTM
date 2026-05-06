# Per-protein Cramér's V analysis
# Teufel et al. (2025) Genome Biology and Evolution — reviewer revision
#
# Computes Cramér's V between broad host category (Human / Avian / Swine)
# and PTM symbol at each site, for each protein.  Results are stratified by
# site rate class (all sites / fast-evolving / slow-evolving).
#
# Outputs written to results/reviewer_analyses/:
#   cramer_v_per_protein.csv  — protein-level mean ± SE
#   cramer_v_per_site.csv     — per-site values with fast/slow flags
#
# Usage (from repo root):
#   python scripts/analysis/cramer_v_by_protein.py

import re, json, csv, warnings
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import chi2_contingency

warnings.filterwarnings("ignore")

N_PERM = 999

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
    # Mammalian_Other
    "Panthera tigris":        "Mammalian_Other",
    "Mus musculus":           "Mammalian_Other",
    "Mustela lutreola":       "Mammalian_Other",
    "Mustela putorius furo":  "Mammalian_Other",
    "Suricata suricatta":     "Mammalian_Other",
    "Felis catus":            "Mammalian_Other",
    # Empty / unknown
    "":                       None,
}


def broad_host(h: str):
    """Map a raw Host string to a broad class using the exact lookup table."""
    return HOST_MAP.get(h.strip() if h else "", None)


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
    """Return list of original alignment positions, one per character column."""
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


def cramers_v(contingency: np.ndarray) -> float:
    if contingency.sum() == 0:
        return np.nan
    chi2 = chi2_contingency(contingency, correction=False)[0]
    n, k = contingency.sum(), min(contingency.shape) - 1
    if k <= 0 or n <= 0:
        return np.nan
    return float(np.sqrt(chi2 / (n * k)))


def chi2_stat(contingency: np.ndarray) -> float:
    """Pearson chi-square statistic for a contingency table."""
    return float(chi2_contingency(contingency, correction=False)[0])


def strain_constrained_perm_p(rows, hosts, symbols, observed_chi2, n_perm=999, rng=None):
    """
    Permutation p-value for host-vs-symbol association while preserving strain structure.

    rows: list of (host, symbol, strain)
    Host labels are permuted only within each strain group.
    """
    if rng is None:
        rng = np.random.default_rng(666)

    h2i = {h: i for i, h in enumerate(hosts)}
    s2i = {s: i for i, s in enumerate(symbols)}

    hosts_arr = np.array([r[0] for r in rows], dtype=object)
    syms_arr = np.array([r[1] for r in rows], dtype=object)
    strains = np.array([r[2] for r in rows], dtype=object)

    by_strain = {st: np.where(strains == st)[0] for st in np.unique(strains)}
    ge = 0

    for _ in range(n_perm):
        perm_hosts = hosts_arr.copy()
        for idx in by_strain.values():
            if len(idx) > 1:
                perm_hosts[idx] = rng.permutation(perm_hosts[idx])

        tbl = np.zeros((len(hosts), len(symbols)), dtype=int)
        for h, s in zip(perm_hosts, syms_arr):
            tbl[h2i[h], s2i[s]] += 1

        perm_chi2 = chi2_stat(tbl)
        if perm_chi2 >= observed_chi2:
            ge += 1

    return float((ge + 1) / (n_perm + 1))


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
meta = load_metadata()
rng = np.random.default_rng(666)
prot_rows, site_rows = [], []

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

    # Resolve taxa to host and strain
    taxon_host = {}
    taxon_strain = {}
    for taxon in records:
        acc  = resolve(taxon)
        info = meta.get(acc, {})
        hb   = info.get("host_broad")
        st   = info.get("strain")
        if hb and st:
            taxon_host[taxon] = hb
            taxon_strain[taxon] = st

    if len(taxon_host) < 10:
        print(f"  skip {protein} — too few resolved taxa ({len(taxon_host)})")
        continue

    sv_all, sv_fast, sv_slow = [], [], []
    n_sig_all = n_sig_fast = n_sig_slow = 0
    n_tested_all = n_tested_fast = n_tested_slow = 0

    for si in range(n_sites):
        char_idx = si + 1
        orig_pos = char_labels[si] if si < len(char_labels) else char_idx

        # Contingency table: host × PTM symbol
        counts = {}
        obs_rows = []
        for taxon, host in taxon_host.items():
            syms = records[taxon]
            if si < len(syms):
                sym = syms[si]
                if sym not in ("-", "?"):
                    counts[(host, sym)] = counts.get((host, sym), 0) + 1
                    obs_rows.append((host, sym, taxon_strain[taxon]))

        if not counts:
            continue

        hosts   = sorted({h for h, s in counts})
        symbols = sorted({s for h, s in counts})
        if len(hosts) < 2 or len(symbols) < 2:
            cv = np.nan
            p_strain_perm = np.nan
            sig_strain_perm = False
        else:
            tbl = np.array([[counts.get((h, s), 0) for s in symbols]
                            for h in hosts])
            cv  = cramers_v(tbl)
            obs_chi2 = chi2_stat(tbl)
            p_strain_perm = strain_constrained_perm_p(
                obs_rows, hosts, symbols, obs_chi2,
                n_perm=N_PERM, rng=rng
            )
            sig_strain_perm = bool(p_strain_perm <= 0.05)

            n_tested_all += 1
            if sig_strain_perm:
                n_sig_all += 1
            if char_idx in fast_idx:
                n_tested_fast += 1
                if sig_strain_perm:
                    n_sig_fast += 1
            if char_idx in slow_idx:
                n_tested_slow += 1
                if sig_strain_perm:
                    n_sig_slow += 1

        site_rows.append({
            "protein":   protein,
            "site_idx":  char_idx,
            "orig_pos":  orig_pos,
            "cramers_v": cv,
            "is_fast":   char_idx in fast_idx,
            "is_slow":   char_idx in slow_idx,
            "n_symbols": len(symbols),
            "n_taxa":    sum(counts.values()),
            "p_strain_perm": p_strain_perm,
            "significant_strain_perm": sig_strain_perm,
        })

        if not np.isnan(cv):
            sv_all.append(cv)
            if char_idx in fast_idx:
                sv_fast.append(cv)
            if char_idx in slow_idx:
                sv_slow.append(cv)

    def mean_se(vals):
        if not vals:
            return np.nan, np.nan
        a = np.array(vals)
        return float(a.mean()), float(a.std() / np.sqrt(len(a)))

    ma, sea = mean_se(sv_all)
    mf, sef = mean_se(sv_fast)
    ms, ses = mean_se(sv_slow)

    prot_rows.append({
        "protein":        protein,
        "n_sites":        n_sites,
        "n_fast_sites":   len(fast_idx),
        "n_slow_sites":   len(slow_idx),
        "n_taxa":         len(taxon_host),
        "cramers_v_all":  round(ma, 4),
        "se_all":         round(sea, 4),
        "cramers_v_fast": round(mf, 4) if not np.isnan(mf) else None,
        "se_fast":        round(sef, 4) if not np.isnan(sef) else None,
        "cramers_v_slow": round(ms, 4) if not np.isnan(ms) else None,
        "se_slow":        round(ses, 4) if not np.isnan(ms) else None,
                "n_sites_tested_strain_perm": n_tested_all,
                "n_sites_sig_strain_perm": n_sig_all,
                "n_fast_tested_strain_perm": n_tested_fast,
                "n_fast_sig_strain_perm": n_sig_fast,
                "n_slow_tested_strain_perm": n_tested_slow,
                "n_slow_sig_strain_perm": n_sig_slow,
    })

    print(f"  {protein:25s}  V_all={ma:.4f}  "
                    f"V_fast={mf:.4f}  V_slow={ms:.4f}  "
                    f"sig_adj={n_sig_all}/{n_tested_all}  "
                    f"(fast_n={len(fast_idx)}, slow_n={len(slow_idx)})")

pd.DataFrame(prot_rows).to_csv(OUT / "cramer_v_per_protein.csv", index=False)
pd.DataFrame(site_rows).to_csv(OUT / "cramer_v_per_site.csv", index=False)
print(f"\nSaved to {OUT}")
