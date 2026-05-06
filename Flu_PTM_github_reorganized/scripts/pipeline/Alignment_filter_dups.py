#!/usr/bin/env python3
"""
Combine per-strain FASTAs, deduplicate, and align with MAFFT
Teufel et al. (2025) Genome Biology and Evolution

Reads one FASTA file per strain directory (H1N1/, H5N1/, H7N9/),
removes sequences with identical content within each strain,
appends the strain name to each sequence header, then runs MAFFT
(--auto) to produce a combined alignment.

Note: the main pipeline (automated_all_6_fixed_alignment_headers.py)
handles deduplication and alignment for all segments automatically.
This standalone script was used for initial single-segment processing.

Usage:
    python Alignment_filter_dups.py [--segment FILENAME] [--mafft PATH]
                                    [--out-combined FILE] [--out-aligned FILE]

Defaults:
    --segment     seg1.fasta.txt
    --mafft       mafft           (assumes mafft is on PATH)
    --out-combined combined_sequences_filtered.fasta
    --out-aligned  aligned_sequences_filtered.fasta
"""

import argparse
import os
import shutil
import subprocess


def parse_fasta(file_path):
    """Read a FASTA file and return {header: sequence} (no gap characters)."""
    sequences = {}
    current_header = None
    current_sequence = []

    with open(file_path) as infile:
        for line in infile:
            line = line.strip()
            if line.startswith(">"):
                if current_header:
                    sequences[current_header] = "".join(current_sequence)
                current_header = line
                current_sequence = []
            elif line:
                current_sequence.append(line)

    if current_header:
        sequences[current_header] = "".join(current_sequence)

    return sequences


def format_sequence(sequence, line_length=60):
    """Wrap a sequence string to the specified line length."""
    return "\n".join(sequence[i:i + line_length] for i in range(0, len(sequence), line_length))


def combine_and_deduplicate(input_dirs, segment_filename, combined_fasta):
    """Combine FASTAs from all strain directories, removing duplicate sequences."""
    with open(combined_fasta, "w") as out:
        for directory in input_dirs:
            filepath = os.path.join(directory, segment_filename)
            if not os.path.exists(filepath):
                print(f"  Skipping {filepath} (not found)")
                continue

            sequences = parse_fasta(filepath)
            seen = {}  # sequence content → header (keeps first occurrence)
            for header, seq in sequences.items():
                if seq not in seen:
                    seen[seq] = header

            for seq, header in seen.items():
                out.write(f"{header}_{directory}\n{format_sequence(seq)}\n")

    print(f"Combined FASTA saved: {combined_fasta}")


def run_mafft(mafft_path, combined_fasta, aligned_fasta):
    """Run MAFFT --auto alignment."""
    cmd = [mafft_path, "--auto", combined_fasta]
    with open(aligned_fasta, "w") as out:
        subprocess.run(cmd, stdout=out, text=True, check=True)
    print(f"Alignment saved: {aligned_fasta}")


def main():
    parser = argparse.ArgumentParser(
        description="Combine per-strain FASTAs, deduplicate, and align with MAFFT."
    )
    parser.add_argument("--segment",      default="seg1.fasta.txt",
                        help="Filename to look for inside each strain directory")
    parser.add_argument("--mafft",        default=shutil.which("mafft") or "mafft",
                        help="Path to the mafft executable (default: mafft on PATH)")
    parser.add_argument("--out-combined", default="combined_sequences_filtered.fasta",
                        help="Output path for the combined deduplicated FASTA")
    parser.add_argument("--out-aligned",  default="aligned_sequences_filtered.fasta",
                        help="Output path for the MAFFT-aligned FASTA")
    args = parser.parse_args()

    input_dirs = ["H1N1", "H5N1", "H7N9"]

    combine_and_deduplicate(input_dirs, args.segment, args.out_combined)
    run_mafft(args.mafft, args.out_combined, args.out_aligned)


if __name__ == "__main__":
    main()