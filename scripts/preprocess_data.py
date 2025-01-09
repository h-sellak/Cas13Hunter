"""
preprocess_data.py

A script for preprocessing raw viral genome sequences for downstream analysis.
This includes removing duplicates, filtering by length, handling ambiguous bases,
calculating GC content, and saving the cleaned sequences.

Author: Hamza SELLAK
Date: 09/01/2025
"""

import os
from collections import Counter
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np


def load_sequences(fasta_file):
    """
    Load viral genome sequences from a FASTA file.
    
    Args:
        fasta_file (str): Path to the FASTA file containing viral genome sequences.
    
    Returns:
        list: A list of SeqRecord objects.
    """
    return list(SeqIO.parse(fasta_file, "fasta"))


def remove_duplicates(sequences):
    """
    Remove duplicate sequences from the dataset.
    
    Args:
        sequences (list): List of SeqRecord objects.
    
    Returns:
        list: A list of unique SeqRecord objects.
    """
    unique_sequences = {str(record.seq): record for record in sequences}
    return list(unique_sequences.values())


def filter_by_length(sequences, min_length, max_length):
    """
    Filter sequences by length.
    
    Args:
        sequences (list): List of SeqRecord objects.
        min_length (int): Minimum acceptable sequence length.
        max_length (int): Maximum acceptable sequence length.
    
    Returns:
        list: A list of SeqRecord objects within the specified length range.
    """
    return [record for record in sequences if min_length <= len(record.seq) <= max_length]


def remove_high_ambiguity(sequences, max_ambiguity_ratio):
    """
    Remove sequences with high ambiguity (excessive 'N' bases).
    
    Args:
        sequences (list): List of SeqRecord objects.
        max_ambiguity_ratio (float): Maximum allowed proportion of ambiguous bases.
    
    Returns:
        list: A list of SeqRecord objects with acceptable ambiguity levels.
    """
    return [
        record for record in sequences
        if str(record.seq).count("N") / len(record.seq) <= max_ambiguity_ratio
    ]


def replace_ambiguous_bases(sequences):
    """
    Replace ambiguous bases ('N') with the most common base at each position.
    
    Args:
        sequences (list): List of SeqRecord objects.
    
    Returns:
        list: A list of SeqRecord objects with ambiguous bases replaced.
    """
    min_length = min(len(record.seq) for record in sequences)
    consensus_seq = "".join(
        Counter([str(record.seq)[i] for record in sequences]).most_common(1)[0][0]
        for i in range(min_length)
    )
    
    def replace_ambiguous(seq, reference):
        seq = list(seq)
        for i, base in enumerate(seq):
            if i < len(reference) and base == "N":
                seq[i] = reference[i]
        return "".join(seq)
    
    return [
        SeqRecord(
            Seq(replace_ambiguous(str(record.seq), consensus_seq)),
            id=record.id,
            name=record.name,
            description=record.description,
            dbxrefs=record.dbxrefs
        )
        for record in sequences
    ]


def calculate_gc_content(seq):
    """
    Calculate GC content of a sequence.
    
    Args:
        seq (str): Nucleotide sequence.
    
    Returns:
        float: GC content as a percentage.
    """
    return (seq.count("G") + seq.count("C")) / len(seq) * 100


def filter_by_gc_content(sequences, min_gc, max_gc):
    """
    Filter sequences by GC content.
    
    Args:
        sequences (list): List of SeqRecord objects.
        min_gc (float): Minimum acceptable GC content percentage.
        max_gc (float): Maximum acceptable GC content percentage.
    
    Returns:
        list: A list of SeqRecord objects within the specified GC content range.
    """
    return [
        record for record in sequences
        if min_gc <= calculate_gc_content(str(record.seq)) <= max_gc
    ]


def save_sequences(sequences, output_file):
    """
    Save processed sequences to a FASTA file.
    
    Args:
        sequences (list): List of SeqRecord objects.
        output_file (str): Path to save the output FASTA file.
    """
    SeqIO.write(sequences, output_file, "fasta")
    print(f"Processed sequences saved to: {output_file}")


if __name__ == "__main__":
    # Input and output file paths
    RAW_DATA_DIR = "../data/raw/"
    PROCESSED_DATA_DIR = "../data/processed/"
    os.makedirs(PROCESSED_DATA_DIR, exist_ok=True)

    fasta_file = os.path.join(RAW_DATA_DIR, "sample_viral_sequences.fasta")
    processed_file = os.path.join(PROCESSED_DATA_DIR, "cleaned_sequences.fasta")

    # Parameters
    MIN_LENGTH = 29000
    MAX_LENGTH = 30000
    MAX_AMBIGUITY_RATIO = 0.05
    MIN_GC = 37
    MAX_GC = 39

    # Processing pipeline
    print("Loading sequences...")
    sequences = load_sequences(fasta_file)
    print(f"Total sequences loaded: {len(sequences)}")

    print("Removing duplicate sequences...")
    sequences = remove_duplicates(sequences)
    print(f"Sequences after removing duplicates: {len(sequences)}")

    print("Filtering sequences by length...")
    sequences = filter_by_length(sequences, MIN_LENGTH, MAX_LENGTH)
    print(f"Sequences after length filtering: {len(sequences)}")

    print("Removing high-ambiguity sequences...")
    sequences = remove_high_ambiguity(sequences, MAX_AMBIGUITY_RATIO)
    print(f"Sequences after removing high-ambiguity sequences: {len(sequences)}")

    print("Replacing ambiguous bases...")
    sequences = replace_ambiguous_bases(sequences)
    print(f"Sequences after replacing ambiguous bases: {len(sequences)}")

    print("Filtering sequences by GC content...")
    sequences = filter_by_gc_content(sequences, MIN_GC, MAX_GC)
    print(f"Sequences after GC content filtering: {len(sequences)}")

    print("Saving processed sequences...")
    save_sequences(sequences, processed_file)
