#!/usr/bin/env python3
"""
Standalone MusiteDeep PTM predictor
Teufel et al. (2025) Genome Biology and Evolution

Queries the MusiteDeep REST API for all 13 PTM types for each sequence
in a FASTA file and writes raw JSON results to an output text file.

Note: the main pipeline (automated_all_6_fixed_alignment_headers.py)
handles PTM prediction as part of the full analysis workflow. This
script provides a lightweight standalone option.

Usage:
    python Get_PMT_from_net.py <input_fasta> <output_txt>

Example:
    python Get_PMT_from_net.py sequences.fasta ptm_predictions.txt
"""

import argparse
import json
import sys
import time

import requests
from Bio import SeqIO

# ---------------------------------------------------------------------------
# MusiteDeep model list (all 13 PTM types used in the paper)
# ---------------------------------------------------------------------------
PTM_MODELS = [
    "Phosphoserine_Phosphothreonine",
    "Phosphotyrosine",
    "N-linked_glycosylation",
    "O-linked_glycosylation",
    "Ubiquitination",
    "SUMOylation",
    "N6-acetyllysine",
    "Methylarginine",
    "Methyllysine",
    "Pyrrolidone_carboxylic_acid",
    "S-palmitoyl_cysteine",
    "Hydroxyproline",
    "Hydroxylysine",
]

API_BASE  = "https://api.musite.net/musitedeep"
API_DELAY = 1  # seconds between requests (be polite to the server)


def predict_ptms(input_fasta, output_txt):
    model_str = ";".join(PTM_MODELS)
    sequences = list(SeqIO.parse(input_fasta, "fasta"))
    print(f"Loaded {len(sequences)} sequences from {input_fasta}")

    with open(output_txt, "w") as out:
        for seq_record in sequences:
            print(f"Processing: {seq_record.id}")

            # Remove alignment gaps before submitting to API
            seq = str(seq_record.seq).replace("-", "")
            url = f"{API_BASE}/{model_str}/{seq}"

            response = requests.get(url)

            if response.ok:
                data = response.json()
                if "Error" in data:
                    print(f"  API error for {seq_record.id}: {data['Error']}")
                else:
                    out.write(f"Sequence: {seq_record.id}\n")
                    out.write(json.dumps(data, indent=4))
                    out.write("\n\n")
            else:
                print(f"  Request failed for {seq_record.id} (HTTP {response.status_code})")
                response.raise_for_status()

            time.sleep(API_DELAY)

    print(f"Done. Results saved to {output_txt}")


def main():
    parser = argparse.ArgumentParser(description="Query MusiteDeep REST API for PTM predictions.")
    parser.add_argument("input_fasta", help="Input FASTA file (aligned or unaligned; gaps are stripped)")
    parser.add_argument("output_txt",  help="Output file for raw JSON results")
    args = parser.parse_args()
    predict_ptms(args.input_fasta, args.output_txt)


if __name__ == "__main__":
    main()