from Bio import AlignIO, SeqIO
from Bio.Align import AlignInfo
from Bio.PDB import PDBParser
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import os
import requests
import glob

def create_consensus_sequence(alignment_file, threshold=0.7):
    """
    Create a consensus sequence from a multiple sequence alignment.
    """
    # Parse the alignment file
    alignment = AlignIO.read(alignment_file, "fasta")
    
    # Create a summary of the alignment
    summary_align = AlignInfo.SummaryInfo(alignment)
    
    # Get the consensus sequence with the specified threshold
    consensus = summary_align.dumb_consensus(threshold=threshold, ambiguous='X')
    
    # Remove any 'X' characters (which represent positions below the threshold)
    # and replace them with the most common amino acid at that position
    consensus_list = list(str(consensus))
    for i, char in enumerate(consensus_list):
        if char == 'X':
            # Get the most frequent amino acid at this position
            column = alignment[:, i]
            aa_counts = {}
            for aa in column:
                if aa != '-':  # Skip gaps
                    aa_counts[aa] = aa_counts.get(aa, 0) + 1
            
            if aa_counts:  # If there are any non-gap amino acids
                most_common_aa = max(aa_counts, key=aa_counts.get)
                consensus_list[i] = most_common_aa
            else:
                # If the column is all gaps, use a gap in the consensus
                consensus_list[i] = '-'
    
    # Convert back to a string
    consensus_str = ''.join(consensus_list)
    
    # Remove gaps for the final sequence
    consensus_str = consensus_str.replace('-', '')
    
    # Create a SeqRecord for the consensus
    consensus_record = SeqRecord(
        Seq(consensus_str),
        id="consensus",
        description="Consensus sequence derived from alignment"
    )
    
    return consensus_record

def predict_structure_esmfold(sequence, output_dir="esmfold_results"):
    """
    Predict the structure of a protein sequence using the ESMFold API.
    """
    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Use the ESM Atlas API for structure prediction
    url = "https://api.esmatlas.com/foldSequence/v1/pdb/"
    
    print(f"Submitting sequence to ESMFold API (length: {len(sequence)})")
    response = requests.post(url, data=sequence)
    
    if response.status_code == 200:
        # Save the predicted structure to a PDB file
        pdb_file = os.path.join(output_dir, "predicted_structure.pdb")
        with open(pdb_file, "w") as f:
            f.write(response.text)
        
        print(f"Structure prediction saved to {pdb_file}")
        return pdb_file
    else:
        print(f"Error: ESMFold API returned status code {response.status_code}")
        print(f"Response: {response.text}")
        return None

def map_alignment_to_consensus_structure(alignment_file, pdb_file, positions_of_interest, output_script="highlight_residues.py"):
    """
    Map positions in an alignment to residues in a consensus structure and create a PyMOL script.
    """
    # Parse the alignment file
    alignment = AlignIO.read(alignment_file, "fasta")
    
    # Create consensus sequence if not already in the alignment
    consensus_record = create_consensus_sequence(alignment_file)
    
    # Get the ungapped consensus sequence
    consensus_seq = str(consensus_record.seq)
    
    # Need to add the consensus sequence to the alignment to map positions
    # Create a gapped version of the consensus to match the alignment
    aligned_consensus = ""
    ungapped_pos = 0
    
    # For each column in the alignment
    for col_idx in range(alignment.get_alignment_length()):
        column = alignment[:, col_idx]
        # Count non-gap characters
        non_gaps = sum(1 for aa in column if aa != '-')
        
        if non_gaps > 0:  # If this column has non-gap characters
            # Add the corresponding consensus character
            aligned_consensus += consensus_seq[ungapped_pos]
            ungapped_pos += 1
        else:
            # Add a gap
            aligned_consensus += '-'
    
    # Create a mapping from alignment positions to consensus sequence positions
    align_pos_to_consensus_pos = {}
    consensus_pos = 1  # 1-indexed positions in the consensus sequence
    
    for align_pos, char in enumerate(aligned_consensus, 1):  # 1-indexed positions in alignment
        if char != '-':  # Not a gap
            align_pos_to_consensus_pos[align_pos] = consensus_pos
            consensus_pos += 1
    
    # Translate alignment positions to consensus sequence positions
    consensus_residues = []
    for pos in positions_of_interest:
        if pos in align_pos_to_consensus_pos:
            consensus_residues.append(align_pos_to_consensus_pos[pos])
        else:
            print(f"Warning: Position {pos} in the alignment does not map to a residue in the consensus sequence.")
    
    # Create a PyMOL script to highlight the residues
    with open(output_script, "w") as f:
        f.write(f"# PyMOL script to highlight residues of interest in the consensus structure\n")
        f.write(f"load {pdb_file}\n")
        f.write(f"hide everything\n")
        f.write(f"show cartoon\n")
        f.write(f"color gray80, all\n")
        
        # Create a selection of residues of interest
        if consensus_residues:
            residue_list = "+".join([str(res) for res in consensus_residues])
            f.write(f"select residues_of_interest, resi {residue_list}\n")
            f.write(f"color red, residues_of_interest\n")
            f.write(f"show sticks, residues_of_interest\n")
            f.write(f"zoom residues_of_interest\n")
            f.write(f"label residues_of_interest and name CA, '%s%s' % (resn, resi)\n")
        else:
            print("No residues of interest were mapped. PyMOL script will not highlight any residues.")
    
    print(f"PyMOL script created: {output_script}")
    print(f"To use it, run: pymol {output_script}")
    
    # Print the mapping for reference
    print("\nAlignment position to consensus residue mapping:")
    for align_pos, cons_pos in sorted(align_pos_to_consensus_pos.items()):
        print(f"Alignment position {align_pos} -> Consensus residue {cons_pos}")
    
    return align_pos_to_consensus_pos

# Main execution code - you can modify these values directly
output_dir = "esmfold_output"
threshold = 0.7
positions_of_interest = [12, 25, 36]  # Edit these values as needed

# Find the first FASTA alignment file in the current directory
alignment_files = glob.glob("*.fasta.txt") + glob.glob("*.fa") + glob.glob("*.aln")

if not alignment_files:
    print("No alignment files found in the current directory.")
else:
    alignment_file = alignment_files[0]
    print(f"Using alignment file: {alignment_file}")
    
    # Create a consensus sequence
    consensus_record = create_consensus_sequence(alignment_file, threshold)
    print(consensus_record)
    # Save the consensus sequence to a file
    consensus_file = os.path.join(output_dir, "consensus.fasta")
    os.makedirs(output_dir, exist_ok=True)
    SeqIO.write(consensus_record, consensus_file, "fasta")
    
    print(f"Consensus sequence saved to {consensus_file}")
    
    # Predict the structure using ESMFold
    pdb_file = predict_structure_esmfold(str(consensus_record.seq), output_dir)
    
    if pdb_file:
        # Map alignment positions to the structure
        mapping = map_alignment_to_consensus_structure(alignment_file, pdb_file, positions_of_interest)