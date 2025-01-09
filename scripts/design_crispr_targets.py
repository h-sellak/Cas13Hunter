"""
design_crispr_targets.py

A script to design CRISPR-Cas13 guide RNAs targeting conserved regions of viral genomes.
This includes:
1. Loading conserved positions and aligned sequences.
2. Extracting sequences around conserved positions.
3. Filtering and evaluating guide RNAs based on length and GC content.
4. Saving valid guide RNAs for downstream analysis.

Author: Hamza SELLAK
Date: 09/01/2025
"""

import os
import numpy as np
from Bio import SeqIO
import pandas as pd


def load_conserved_positions(file_path):
    """
    Load conserved positions from a NumPy file.

    Args:
        file_path (str): Path to the .npy file containing conserved positions.

    Returns:
        numpy.ndarray: Array of conserved positions.
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"Conserved positions file not found: {file_path}")
    return np.load(file_path)


def load_aligned_sequences(file_path):
    """
    Load aligned sequences from a FASTA file.

    Args:
        file_path (str): Path to the FASTA file containing aligned sequences.

    Returns:
        list: List of SeqRecord objects from the alignment.
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"Aligned sequences file not found: {file_path}")
    return list(SeqIO.parse(file_path, "fasta"))


def extract_guide_candidates(conserved_positions, alignment, window_size):
    """
    Extract sequences around conserved positions as guide RNA candidates.

    Args:
        conserved_positions (numpy.ndarray): Array of conserved positions.
        alignment (list): List of SeqRecord objects from the alignment.
        window_size (int): Number of bases around conserved positions to extract.

    Returns:
        list: List of guide RNA candidates with metadata.
    """
    guide_candidates = []
    for pos in conserved_positions:
        for record in alignment:
            start = max(0, pos - window_size // 2)
            end = min(len(record.seq), pos + window_size // 2)
            guide_seq = str(record.seq[start:end]).replace("-", "")  # Remove gaps
            if len(guide_seq) == window_size:
                guide_candidates.append({
                    "conserved_pos": pos,
                    "sequence_id": record.id,
                    "guide_sequence": guide_seq
                })
    return guide_candidates


def calculate_gc_content(seq):
    """
    Calculate GC content of a nucleotide sequence.

    Args:
        seq (str): Nucleotide sequence.

    Returns:
        float: GC content as a percentage.
    """
    seq = seq.upper()
    return (seq.count("G") + seq.count("C")) / len(seq) * 100 if len(seq) > 0 else 0


def evaluate_guides(guide_candidates, guide_length, gc_range):
    """
    Filter and evaluate guide RNA candidates based on length and GC content.

    Args:
        guide_candidates (list): List of guide RNA candidates.
        guide_length (int): Desired guide RNA length.
        gc_range (tuple): Acceptable GC content range (min_gc, max_gc).

    Returns:
        list: List of valid guide RNAs with metadata.
    """
    valid_guides = []
    for guide in guide_candidates:
        guide_seq = guide["guide_sequence"][:guide_length].upper()
        if len(guide_seq) == guide_length:
            gc_content = calculate_gc_content(guide_seq)
            if gc_range[0] <= gc_content <= gc_range[1]:
                valid_guides.append({
                    "conserved_pos": guide["conserved_pos"],
                    "sequence_id": guide["sequence_id"],
                    "guide_sequence": guide_seq,
                    "gc_content": gc_content
                })
    return valid_guides


def save_guides(guides, output_file):
    """
    Save valid guide RNAs to a CSV file.

    Args:
        guides (list): List of valid guide RNAs with metadata.
        output_file (str): Path to the output CSV file.

    Returns:
        None
    """
    df = pd.DataFrame(guides)
    df.to_csv(output_file, index=False)
    print(f"Designed guides saved to: {output_file}")


if __name__ == "__main__":
    # Define file paths
    PROCESSED_DATA_DIR = "../data/processed/"
    conserved_file = os.path.join(PROCESSED_DATA_DIR, "conserved_regions.npy")
    aligned_file = os.path.join(PROCESSED_DATA_DIR, "aligned_sequences.fasta")
    guides_file = os.path.join(PROCESSED_DATA_DIR, "designed_guides.csv")

    # Parameters
    WINDOW_SIZE = 30  # Number of bases around conserved positions
    GUIDE_LENGTH = 22  # Desired guide RNA length
    GC_RANGE = (30, 70)  # Acceptable GC content range in percentage

    # Load conserved positions and aligned sequences
    print("Loading conserved positions...")
    conserved_positions = load_conserved_positions(conserved_file)
    print(f"Loaded {len(conserved_positions)} conserved positions.")

    print("Loading aligned sequences...")
    alignment = load_aligned_sequences(aligned_file)
    print(f"Loaded {len(alignment)} aligned sequences.")

    # Extract guide RNA candidates
    print("Extracting guide RNA candidates...")
    guide_candidates = extract_guide_candidates(conserved_positions, alignment, WINDOW_SIZE)
    print(f"Extracted {len(guide_candidates)} guide RNA candidates.")

    # Evaluate guide RNA candidates
    print("Evaluating guide RNA candidates...")
    valid_guides = evaluate_guides(guide_candidates, GUIDE_LENGTH, GC_RANGE)
    print(f"Generated {len(valid_guides)} valid guide RNAs.")

    # Save valid guides
    print("Saving valid guide RNAs...")
    save_guides(valid_guides, guides_file)
