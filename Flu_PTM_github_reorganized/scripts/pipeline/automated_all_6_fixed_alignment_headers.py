#!/usr/bin/env python3
"""
IAV PTM analysis pipeline
Teufel et al. (2025) Genome Biology and Evolution

For each of the 8 IAV genomic segments across H1N1, H5N1, and H7N9,
this pipeline:
  1. Extracts protein sequences from per-strain FASTA files
  2. Aligns sequences with MAFFT (--auto)
  3. Builds a maximum-likelihood tree with IQ-TREE 2 (ModelFinder Plus,
     1000 ultrafast bootstraps)
  4. Predicts PTM sites via the MusiteDeep REST API (13 PTM types,
     confidence threshold 0.75)
  5. Constructs the discrete character matrix (*_multistate.txt) for RevBayes
  6. Writes position-mapping tables (alignment ↔ protein coordinates)

Usage:
    python automated_all_6_fixed_alignment_headers.py [--output-dir DIR]
                                                       [--mafft PATH]
                                                       [--iqtree PATH]
                                                       [--force-realign]

Input layout expected:
    <strain>/<segment>.fasta.txt   for strain in {H1N1, H5N1, H7N9}
                                   for segment in {seg1 … seg8}
"""

import argparse
import json
import logging
import os
import re
import subprocess
import sys
import time
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed

import requests
from Bio import AlignIO, SeqIO


# ---------------------------------------------------------------------------
# PTM model string (MusiteDeep API parameter)
# ---------------------------------------------------------------------------
PTM_MODELS = ";".join([
    "Phosphoserine_Phosphothreonine", "Phosphotyrosine",
    "N-linked_glycosylation", "O-linked_glycosylation",
    "Ubiquitination", "SUMOylation", "N6-acetyllysine",
    "Methylarginine", "Methyllysine",
    "Pyrrolidone_carboxylic_acid", "S-palmitoyl_cysteine",
    "Hydroxyproline", "Hydroxylysine",
])

PTM_CONFIDENCE_THRESHOLD = 0.75

# Map PTM type names to single-character symbols used in the multistate matrix
PTM_SYMBOLS = {
    "Phosphoserine_Phosphothreonine": "S",
    "Phosphotyrosine":                "Y",
    "N-linked_glycosylation":         "N",
    "O-linked_glycosylation":         "O",
    "Ubiquitination":                 "K",
    "SUMOylation":                    "U",
    "N6-acetyllysine":                "A",
    "Methylarginine":                 "R",
    "Methyllysine":                   "M",
    "Pyrrolidone_carboxylic_acid":    "Q",
    "S-palmitoyl_cysteine":           "C",
    "Hydroxyproline":                 "H",
    "Hydroxylysine":                  "L",
}

UNMODIFIED_SYMBOL = "X"
GAP_SYMBOL        = "-"

# ---------------------------------------------------------------------------
# Segment → protein name mapping
# Keys are substrings that appear in FASTA description lines.
# ---------------------------------------------------------------------------
SEGMENT_PROTEIN_PATTERNS = {
    "seg1": {
        "|PB2|": "PB2_polymerase",
        "|polymerase PB2|": "PB2_polymerase",
        "|polymerase basic protein 2|": "PB2_polymerase",
    },
    "seg2": {
        "|PB1|": "PB1_polymerase",
        "|polymerase PB1|": "PB1_polymerase",
        "|polymerase basic protein 1|": "PB1_polymerase",
        "|PB1-F2|": "PB1-F2_protein",
        "|polymerase basic protein 1-F2|": "PB1-F2_protein",
        "|PB1-N40|": "PB1-N40_protein",
    },
    "seg3": {
        "|PA|": "PA_polymerase",
        "|polymerase PA|": "PA_polymerase",
        "|polymerase acidic protein|": "PA_polymerase",
        "|PA-X|": "PA-X_protein",
        "|polymerase PA-X protein|": "PA-X_protein",
    },
    "seg4": {
        "|hemagglutinin|": "HA_full",
        "|HA|": "HA_full",
        "|HA1|": "HA1_domain",
        "|HA2|": "HA2_domain",
    },
    "seg5": {
        "|NP|": "NP_protein",
        "|nucleoprotein|": "NP_protein",
        "|nucleocapsid protein|": "NP_protein",
    },
    "seg6": {
        "|NA|": "NA_neuraminidase",
        "|neuraminidase|": "NA_neuraminidase",
    },
    "seg7": {
        "|M1|": "M1_matrix_protein",
        "|matrix protein 1|": "M1_matrix_protein",
        "|M2|": "M2_ion_channel",
        "|matrix protein 2|": "M2_ion_channel",
        "|ion channel protein M2|": "M2_ion_channel",
    },
    "seg8": {
        "|NS1|": "NS1_protein",
        "|nonstructural protein 1|": "NS1_protein",
        "|NS2|": "NS2_protein",
        "|NEP|": "NEP_protein",
        "|nuclear export protein|": "NEP_protein",
        "|nonstructural protein 2|": "NS2_protein",
    },
}


# ---------------------------------------------------------------------------
# Pipeline class
# ---------------------------------------------------------------------------

class IAVPipeline:
    """End-to-end pipeline from raw FASTA to RevBayes input."""

    def __init__(self, strains, output_dir, mafft_path, iqtree_path,
                 force_realign=False):
        self.strains       = strains
        self.output_dir    = output_dir
        self.mafft_path    = mafft_path
        self.iqtree_path   = iqtree_path
        self.force_realign = force_realign

        logging.basicConfig(
            level=logging.INFO,
            format="%(asctime)s %(levelname)s %(message)s",
            handlers=[
                logging.FileHandler("pipeline.log"),
                logging.StreamHandler(sys.stdout),
            ],
        )
        self.log = logging.getLogger(__name__)

    # ------------------------------------------------------------------
    # Step 1: Extract protein sequences from segment FASTA files
    # ------------------------------------------------------------------

    def extract_proteins(self, segment):
        """Return {protein_name: [(header, sequence), ...]} for one segment."""
        proteins = defaultdict(list)
        patterns = SEGMENT_PROTEIN_PATTERNS.get(segment, {})

        for strain in self.strains:
            fasta_path = os.path.join(strain, f"{segment}.fasta.txt")
            if not os.path.exists(fasta_path):
                self.log.warning("Missing: %s", fasta_path)
                continue

            for record in SeqIO.parse(fasta_path, "fasta"):
                desc = record.description
                for pattern, protein_name in patterns.items():
                    if pattern in desc:
                        safe_id = re.sub(r"[^\w\-|]", "_", record.id)
                        proteins[protein_name].append(
                            (f"{strain}__{safe_id}", str(record.seq))
                        )
                        break

        return proteins

    # ------------------------------------------------------------------
    # Step 2: Multiple sequence alignment (MAFFT)
    # ------------------------------------------------------------------

    def align(self, protein_name, sequences, segment_dir):
        """Align sequences with MAFFT --auto. Returns path to aligned FASTA."""
        raw_fasta     = os.path.join(segment_dir, f"{protein_name}_sequences.fasta")
        aligned_fasta = os.path.join(segment_dir, f"{protein_name}_aligned.fasta")

        if os.path.exists(aligned_fasta) and not self.force_realign:
            self.log.info("Using cached alignment: %s", aligned_fasta)
            return aligned_fasta

        with open(raw_fasta, "w") as fh:
            for header, seq in sequences:
                fh.write(f">{header}\n{seq}\n")

        result = subprocess.run(
            [self.mafft_path, "--auto", raw_fasta],
            capture_output=True, text=True
        )
        if result.returncode != 0:
            self.log.error("MAFFT failed:\n%s", result.stderr)
            return None

        with open(aligned_fasta, "w") as fh:
            fh.write(result.stdout)

        self.log.info("Alignment complete: %s", aligned_fasta)
        return aligned_fasta

    # ------------------------------------------------------------------
    # Step 3: Maximum-likelihood tree (IQ-TREE 2)
    # ------------------------------------------------------------------

    def build_tree(self, aligned_fasta):
        """Run IQ-TREE 2. Returns path to NEXUS tree file."""
        tree_nex = aligned_fasta.replace("_aligned.fasta", "_tree.nex")

        if os.path.exists(tree_nex):
            self.log.info("Using cached tree: %s", tree_nex)
            return tree_nex

        cmd = [
            self.iqtree_path, "-s", aligned_fasta,
            "-m", "MFP", "-bb", "1000", "-nt", "AUTO",
            "--prefix", aligned_fasta.replace("_aligned.fasta", "_iqtree"),
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            self.log.error("IQ-TREE failed:\n%s", result.stderr)
            return None

        treefile = aligned_fasta.replace("_aligned.fasta", "_iqtree.treefile")
        if not os.path.exists(treefile):
            self.log.error("IQ-TREE output not found: %s", treefile)
            return None

        # Build tip-name mapping (full FASTA ID → accession)
        id_map = {}
        for record in SeqIO.parse(aligned_fasta, "fasta"):
            m = re.search(r"([A-Z]{1,2}_?\d{6,})", record.id)
            id_map[record.id] = m.group(1) if m else record.id[:20]

        mapping_file = aligned_fasta.replace("_aligned.fasta", "_id_mapping.json")
        with open(mapping_file, "w") as fh:
            json.dump(id_map, fh, indent=2)

        # Convert to NEXUS with shortened tip labels
        with open(treefile) as fh:
            newick = fh.read().strip()
        for long_id, short_id in id_map.items():
            newick = newick.replace(long_id, short_id)

        n_taxa = len(id_map)
        nexus = (
            "#NEXUS\n\nBEGIN TAXA;\n"
            f"\tDIMENSIONS ntax={n_taxa};\n"
            f"\tTAXLABELS {' '.join(id_map.values())};\n"
            "END;\n\nBEGIN TREES;\n"
            f"\tTREE tree_1 = {newick}\n"
            "END;\n"
        )
        with open(tree_nex, "w") as fh:
            fh.write(nexus)

        self.log.info("NEXUS tree saved: %s", tree_nex)
        return tree_nex

    # ------------------------------------------------------------------
    # Step 4: PTM prediction via MusiteDeep REST API
    # ------------------------------------------------------------------

    def _position_map(self, aligned_seq):
        """Return aligned→unaligned and unaligned→aligned index dicts."""
        aln2unaln, unaln2aln = {}, {}
        u = 0
        for a, aa in enumerate(aligned_seq):
            if aa != "-":
                aln2unaln[a] = u
                unaln2aln[u] = a
                u += 1
        return aln2unaln, unaln2aln

    def _query_musite(self, short_id, ungapped_seq, retries=3):
        """Query MusiteDeep API with exponential back-off. Returns JSON or None."""
        url   = f"https://api.musite.net/musitedeep/{PTM_MODELS}/{ungapped_seq}"
        delay = 5
        for attempt in range(retries):
            try:
                r = requests.get(url, timeout=60)
                if r.ok:
                    return r.json()
                self.log.warning("HTTP %d for %s (attempt %d/%d)",
                                 r.status_code, short_id, attempt + 1, retries)
            except Exception as exc:
                self.log.warning("Request error for %s: %s (attempt %d/%d)",
                                 short_id, exc, attempt + 1, retries)
            time.sleep(delay)
            delay *= 2
        return None

    def predict_ptms(self, aligned_fasta, id_map):
        """
        Predict PTMs for all sequences. Returns a dict:
            {short_id: {aligned_position: ptm_symbol}}
        and writes a raw predictions text file.
        """
        ptm_file = aligned_fasta.replace("_aligned.fasta", "_ptm_predictions.txt")
        cache    = aligned_fasta.replace("_aligned.fasta", "_alignment_data.json")

        if os.path.exists(cache) and not self.force_realign:
            self.log.info("Using cached PTM data: %s", cache)
            with open(cache) as fh:
                return json.load(fh)

        records = list(SeqIO.parse(aligned_fasta, "fasta"))

        def _worker(record):
            short_id   = id_map.get(record.id, record.id)
            ungapped   = str(record.seq).replace("-", "")
            aln2u, u2a = self._position_map(str(record.seq))
            if len(ungapped) < 10:
                return short_id, {}, aln2u
            data = self._query_musite(short_id, ungapped)
            return short_id, data or {}, aln2u

        results = {}
        with ThreadPoolExecutor(max_workers=4) as pool:
            futures = {pool.submit(_worker, r): r.id for r in records}
            for fut in as_completed(futures):
                short_id, data, aln2u = fut.result()
                results[short_id] = {"ptm_data": data, "aln_to_ungapped": aln2u}

        # Write raw predictions
        with open(ptm_file, "w") as fh:
            for sid, entry in results.items():
                fh.write(f">{sid}\n{json.dumps(entry['ptm_data'])}\n")

        with open(cache, "w") as fh:
            json.dump(results, fh)

        self.log.info("PTM predictions saved: %s", ptm_file)
        return results

    # ------------------------------------------------------------------
    # Step 5: Build RevBayes multistate character matrix
    # ------------------------------------------------------------------

    @staticmethod
    def _ptm_symbol_at(position_ungapped, ptm_data):
        """
        Return the highest-confidence PTM symbol at a given unaligned
        position, or UNMODIFIED_SYMBOL if nothing exceeds the threshold.
        """
        if not ptm_data:
            return UNMODIFIED_SYMBOL

        # MusiteDeep returns a list of {position, type, score, ...} dicts
        hits = []
        for entry in ptm_data:
            pos   = entry.get("position", -1) - 1  # convert to 0-based
            score = entry.get("score", 0.0)
            ptype = entry.get("type", "")
            if pos == position_ungapped and score >= PTM_CONFIDENCE_THRESHOLD:
                hits.append((score, PTM_SYMBOLS.get(ptype, UNMODIFIED_SYMBOL)))
        if not hits:
            return UNMODIFIED_SYMBOL
        return max(hits, key=lambda x: x[0])[1]

    def build_multistate_matrix(self, aligned_fasta, ptm_results, id_map):
        """
        Construct the NEXUS-format discrete character matrix for RevBayes.
        Each column is an alignment position; the state at a position is
        the PTM symbol (or X = unmodified, - = gap).
        """
        out_file = aligned_fasta.replace("_aligned.fasta", "_multistate.txt")

        records    = list(SeqIO.parse(aligned_fasta, "fasta"))
        aln_length = len(records[0].seq)
        n_seqs     = len(records)

        rows = []
        for record in records:
            short_id   = id_map.get(record.id, record.id)
            aln_seq    = str(record.seq)
            entry      = ptm_results.get(short_id, {})
            aln2u      = entry.get("aln_to_ungapped", {})
            ptm_data   = (entry.get("ptm_data") or {})

            # Flatten MusiteDeep nested structure to a list of prediction dicts
            flat_predictions = []
            if isinstance(ptm_data, dict):
                for ptype, preds in ptm_data.items():
                    if isinstance(preds, list):
                        for p in preds:
                            p["type"] = ptype
                            flat_predictions.append(p)
            elif isinstance(ptm_data, list):
                flat_predictions = ptm_data

            row_chars = []
            for aln_pos in range(aln_length):
                aa = aln_seq[aln_pos]
                if aa == "-":
                    row_chars.append(GAP_SYMBOL)
                else:
                    unaln_pos = aln2u.get(aln_pos)
                    symbol = (
                        self._ptm_symbol_at(unaln_pos, flat_predictions)
                        if unaln_pos is not None else UNMODIFIED_SYMBOL
                    )
                    row_chars.append(symbol)

            rows.append((short_id, "".join(row_chars)))

        all_symbols = sorted(set(c for _, seq in rows for c in seq))
        symbol_str  = " ".join(all_symbols)

        with open(out_file, "w") as fh:
            fh.write("#NEXUS\n\nBEGIN DATA;\n")
            fh.write(f"\tDIMENSIONS ntax={n_seqs} nchar={aln_length};\n")
            fh.write(f'\tFORMAT DATATYPE=STANDARD SYMBOLS="{symbol_str}" '
                     f'MISSING=? GAP={GAP_SYMBOL};\n')
            fh.write("\tMATRIX\n")
            for taxon, seq in rows:
                fh.write(f"\t{taxon:<30} {seq}\n")
            fh.write("\t;\nEND;\n")

        self.log.info("Multistate matrix saved: %s", out_file)
        return out_file

    # ------------------------------------------------------------------
    # Step 6: Write position metadata
    # ------------------------------------------------------------------

    def write_position_metadata(self, aligned_fasta, ptm_results, id_map):
        """Write a tab-delimited table of alignment-position metadata."""
        meta_file = aligned_fasta.replace("_aligned.fasta",
                                           "_position_metadata.txt")
        records    = list(SeqIO.parse(aligned_fasta, "fasta"))
        aln_length = len(records[0].seq)

        # Use the first sequence's ungapped consensus as column reference
        first_id  = id_map.get(records[0].id, records[0].id)
        entry     = ptm_results.get(first_id, {})
        aln2u     = entry.get("aln_to_ungapped", {})

        with open(meta_file, "w") as fh:
            fh.write("aln_position\tungapped_position\n")
            for aln_pos in range(aln_length):
                u_pos = aln2u.get(aln_pos, "")
                fh.write(f"{aln_pos + 1}\t{u_pos + 1 if u_pos != '' else 'gap'}\n")

        self.log.info("Position metadata saved: %s", meta_file)

    # ------------------------------------------------------------------
    # Orchestration: one segment
    # ------------------------------------------------------------------

    def process_segment(self, segment):
        proteins = self.extract_proteins(segment)
        if not proteins:
            self.log.warning("No proteins extracted for %s", segment)
            return

        segment_dir = os.path.join(self.output_dir, segment)
        os.makedirs(segment_dir, exist_ok=True)

        for protein_name, sequences in proteins.items():
            self.log.info("Processing %s / %s (%d sequences)",
                          segment, protein_name, len(sequences))

            aligned = self.align(protein_name, sequences, segment_dir)
            if aligned is None:
                continue

            tree = self.build_tree(aligned)

            # Load the tip-name mapping written by build_tree
            mapping_path = aligned.replace("_aligned.fasta", "_id_mapping.json")
            id_map = {}
            if os.path.exists(mapping_path):
                with open(mapping_path) as fh:
                    id_map = json.load(fh)

            ptm_results = self.predict_ptms(aligned, id_map)
            self.build_multistate_matrix(aligned, ptm_results, id_map)
            self.write_position_metadata(aligned, ptm_results, id_map)

    def run(self):
        for seg in sorted(SEGMENT_PROTEIN_PATTERNS):
            self.log.info("=== Segment %s ===", seg)
            self.process_segment(seg)
        self.log.info("Pipeline complete.")


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def _find_on_path(name):
    """Return the first executable named `name` found on PATH, or None."""
    for directory in os.environ.get("PATH", "").split(os.pathsep):
        candidate = os.path.join(directory, name)
        if os.path.isfile(candidate) and os.access(candidate, os.X_OK):
            return candidate
    return None


def main():
    parser = argparse.ArgumentParser(
        description="IAV PTM pipeline (MAFFT → IQ-TREE → MusiteDeep → RevBayes input)"
    )
    parser.add_argument("--output-dir",   default="processed_segments_fixed_final",
                        help="Directory for all output files")
    parser.add_argument("--mafft",        default=None,
                        help="Path to MAFFT executable (default: auto-detect on PATH)")
    parser.add_argument("--iqtree",       default=None,
                        help="Path to IQ-TREE 2 executable (default: auto-detect on PATH)")
    parser.add_argument("--force-realign", action="store_true",
                        help="Rerun alignment and PTM prediction even if cached files exist")
    args = parser.parse_args()

    mafft  = args.mafft  or _find_on_path("mafft")  or _find_on_path("mafft.bat")
    iqtree = args.iqtree or _find_on_path("iqtree2") or _find_on_path("iqtree2.exe")

    if mafft is None:
        sys.exit("MAFFT not found. Install it or pass --mafft PATH.")
    if iqtree is None:
        sys.exit("IQ-TREE 2 not found. Install it or pass --iqtree PATH.")

    pipeline = IAVPipeline(
        strains       = ["H1N1", "H5N1", "H7N9"],
        output_dir    = args.output_dir,
        mafft_path    = mafft,
        iqtree_path   = iqtree,
        force_realign = args.force_realign,
    )
    pipeline.run()


if __name__ == "__main__":
    main()
