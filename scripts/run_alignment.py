"""
run_alignment.py

A script for performing multiple sequence alignment (MSA) on processed viral genome sequences
using MAFFT with the FFT-NS-2 algorithm for speed. The script automatically determines the 
optimal number of CPU threads to use based on the system's available cores.

Author: Hamza SELLAK
Date: 09/01/2025
"""

import os
from Bio.Align.Applications import MafftCommandline


def get_optimal_threads():
    """
    Determine the optimal number of CPU threads to use for alignment.
    
    Returns:
        int: The number of threads to use, calculated as 75% of available CPU cores.
    """
    num_cores = os.cpu_count()
    return max(1, int(num_cores * 0.75))


def run_mafft(input_file, output_file, num_threads):
    """
    Perform multiple sequence alignment using MAFFT with FFT-NS-2 mode.
    
    Args:
        input_file (str): Path to the input FASTA file containing sequences.
        output_file (str): Path to save the aligned sequences in FASTA format.
        num_threads (int): Number of CPU threads to use for MAFFT.
    
    Returns:
        None
    """
    print(f"Running MAFFT alignment on {input_file} with {num_threads} threads...")
    
    # Set up MAFFT command
    mafft_cline = MafftCommandline(
        input=input_file,
        thread=num_threads,  # Use optimal number of threads
        auto=True            # Automatically chooses FFT-NS-2 for speed
    )
    
    # Execute MAFFT and capture the output
    stdout, stderr = mafft_cline()
    
    # Save the aligned sequences to the output file
    with open(output_file, "w") as f:
        f.write(stdout)
    
    print(f"Alignment complete. Results saved to: {output_file}")


if __name__ == "__main__":
    # Define input and output file paths
    PROCESSED_DATA_DIR = "../data/processed/"
    input_file = os.path.join(PROCESSED_DATA_DIR, "cleaned_sequences.fasta")
    output_file = os.path.join(PROCESSED_DATA_DIR, "aligned_sequences.fasta")
    
    # Ensure the input file exists
    if not os.path.exists(input_file):
        raise FileNotFoundError(f"Input file not found: {input_file}")
    
    # Get the optimal number of threads
    num_threads = get_optimal_threads()
    print(f"Optimal number of threads determined: {num_threads}")
    
    # Run the MAFFT alignment
    run_mafft(input_file, output_file, num_threads)
