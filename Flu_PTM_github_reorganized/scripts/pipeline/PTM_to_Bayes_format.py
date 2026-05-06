#!/usr/bin/env python3
"""
Convert MusiteDeep PTM predictions to RevBayes multistate character matrix
Teufel et al. (2025) Genome Biology and Evolution

Reads all *_ptm_predictions.txt files in a directory (produced by
Get_PMT_from_net.py), selects the highest-scoring PTM at each position,
applies a confidence filter, and writes a whitespace-delimited multistate
matrix (*_multistate.txt) where each column is a PTM-site position and
each row is a sequence. Positions without a PTM above the threshold are
encoded as 'X'; unrecognized PTM types are encoded as '?'.

Note: this is an early standalone script used for initial single-segment
testing. The main pipeline (automated_all_6_fixed_alignment_headers.py)
supersedes it and uses a revised PTM-to-symbol mapping and a stricter
confidence threshold (0.75). The symbol mapping here reflects an earlier
convention and is retained for archival purposes.

Usage:
    python PTM_to_Bayes_format.py [--input-dir DIR] [--min-score FLOAT]

Defaults:
    --input-dir   Seg_2_alignment
    --min-score   0.5
"""

import argparse
import os


# ---------------------------------------------------------------------------
# PTM type → single-character symbol
# Note: this mapping predates the final paper convention used in
# automated_all_6_fixed_alignment_headers.py. In the final analysis,
# Ubiquitination → "K" and SUMOylation → "U".
# ---------------------------------------------------------------------------
PTM_TO_SYMBOL = {
    "Phosphoserine":              "S",
    "Phosphothreonine":           "T",
    "Phosphotyrosine":            "Y",
    "N-linked_glycosylation":     "N",
    "O-linked_glycosylation":     "O",
    "Ubiquitination":             "U",
    "SUMOylation":                "M",
    "N6-acetyllysine":            "A",
    "Methylarginine":             "R",
    "Methyllysine":               "K",
    "Pyrrolidone_carboxylic_acid":"Q",
    "S-palmitoyl_cysteine":       "C",
    "Hydroxyproline":             "H",
    "Hydroxylysine":              "L",
}

UNMODIFIED_SYMBOL = "X"


def parse_ptm_score(ptm_string):
    """Return (ptm_type, score) for the highest-scoring PTM in a semicolon-delimited string."""
    best_score = 0.0
    best_ptm   = None
    for pair in ptm_string.split(";"):
        ptm, score = pair.split(":")
        score = float(score)
        if score > best_score:
            best_score = score
            best_ptm   = ptm
    return best_ptm, best_score


def process_predictions(input_dir, min_score):
    """Convert all *_ptm_predictions.txt files in input_dir to multistate matrices."""
    for filename in os.listdir(input_dir):
        if not filename.endswith("_ptm_predictions.txt"):
            continue

        input_path  = os.path.join(input_dir, filename)
        output_path = os.path.join(input_dir, filename.replace("_ptm_predictions", "_multistate.txt"))

        high_confidence_ptms = {}
        current_protein      = None

        with open(input_path) as f:
            for line in f:
                line = line.strip()

                if line.startswith("Sequence:"):
                    current_protein = line.split("Sequence:", 1)[1].strip()
                    high_confidence_ptms[current_protein] = {}

                elif line.startswith("Raw prediction data:"):
                    if current_protein is None:
                        continue

                    try:
                        # Parse the "fake dict" format written by Get_PMT_from_net.py
                        raw   = line.replace("Raw prediction data:", "").strip().strip("{}")
                        parts = raw.split(", ")
                        data  = {}
                        for part in parts:
                            key, value = part.split(": ", 1)
                            data[key.strip("'\" ")] = value.strip("'\" ")

                        ptm_type, score = parse_ptm_score(data["PTMscores"])
                        position        = int(data["Position"])

                        if score > min_score:
                            high_confidence_ptms[current_protein][position] = ptm_type

                    except Exception as exc:
                        print(f"Warning: could not parse line:\n  {line}\n  {exc}")

        all_positions = sorted(
            {pos for ptms in high_confidence_ptms.values() for pos in ptms}
        )

        with open(output_path, "w") as out:
            for protein_id, ptms in high_confidence_ptms.items():
                out.write(protein_id)
                for pos in all_positions:
                    symbol = PTM_TO_SYMBOL.get(ptms[pos], "?") if pos in ptms else UNMODIFIED_SYMBOL
                    out.write(f" {symbol}")
                out.write("\n")

        print(f"Processed {filename} -> {output_path}")
        print(f"  Positions: {all_positions}")


def main():
    parser = argparse.ArgumentParser(
        description="Convert MusiteDeep prediction files to RevBayes multistate matrices."
    )
    parser.add_argument("--input-dir",  default="Seg_2_alignment",
                        help="Directory containing *_ptm_predictions.txt files")
    parser.add_argument("--min-score",  type=float, default=0.5,
                        help="Minimum confidence score to retain a PTM call (default: 0.5)")
    args = parser.parse_args()
    process_predictions(args.input_dir, args.min_score)


if __name__ == "__main__":
    main()
