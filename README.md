# Cas13Hunter

A bioinformatics pipeline for identifying conserved RNA regions in viral genomes and designing CRISPR-Cas13 guide RNAs for therapeutic development. This tool focuses on enabling rapid target identification for RNA respiratory viruses like SARS-CoV-2, leveraging multiple sequence alignment, functional annotation, and target prediction techniques.

## Getting Started

These instructions will help you set up the project environment and get started with Cas13Hunter.

---

## Prerequisites

To run this project, ensure you have the following installed:

1. [Anaconda/Miniconda](https://docs.conda.io/en/latest/miniconda.html)
2. [Git](https://git-scm.com/)
3. [MAFFT](https://mafft.cbrc.jp/alignment/software/)

---

## Setting Up the Environment

1. **Clone the Repository**:
   ```bash
   git clone https://github.com/h-sellak/Cas13Hunter.git
   cd Cas13Hunter
   ```

2. **Create and Activate the Conda Environment**:
   - Create the environment from the provided `environment.yml` file:
     ```bash
     conda env create -f environment.yml
     ```
   - Activate the environment:
     ```bash
     conda activate cas13hunter
     ```

3. **Verify Installation**:
   Ensure all required packages are installed by listing them:
   ```bash
   conda list
   ```

---

## Directory Structure

```
Cas13Hunter/
├── data/
│   ├── raw/               # Raw viral genome sequences
│   ├── processed/         # Processed and aligned sequences
├── notebooks/             # Jupyter notebooks for experimentation
│   ├── 01_data_preprocessing.ipynb     # Preprocessing of raw sequences
│   ├── 02_alignment_analysis.ipynb     # Multiple sequence alignment and conserved region analysis
│   ├── 03_crispr_target_design.ipynb   # CRISPR-Cas13 guide RNA design
├── scripts/               # Python scripts for pipeline modules
│   ├── preprocess_data.py             # Script for sequence preprocessing
│   ├── run_alignment.py               # Script for running MAFFT alignment
│   ├── design_crispr_targets.py       # Script to design CRISPR-Cas13 guide RNAs
├── environment.yml        # Conda environment configuration file
└── README.md              # Project documentation

```

---

## Notebooks

1. **`01_data_preprocessing.ipynb`:**

   - Prepares raw viral genome sequences for analysis by:
      - Removing duplicates.
      - Filtering sequences by length and ambiguity.
      - Calculating and filtering by GC content.
      - Replacing ambiguous bases with likely alternatives.
   - Output: Cleaned sequences saved to `data/processed/cleaned_sequences.fasta`.

2. **`02_alignment_analysis.ipynb`:**

   - Performs multiple sequence alignment (MSA) using MAFFT.
   - Analyzes conservation scores to identify conserved genome regions.
   - Visualizes conservation across genome positions.
   - Output: Aligned sequences saved to `data/processed/aligned_sequences.fasta` and conserved positions saved to `data/processed/conserved_regions.npy`.


---

## Scripts

`preprocess_data.py`
   
   - Preprocesses raw viral genome sequences by:
      - Removing duplicates.
      - Filtering by length and ambiguity.
      - Replacing ambiguous bases.
      - Filtering sequences by GC content.
   - Usage:

   ```bash
   python scripts/preprocess_data.py
   ```

`run_alignment.py`
   
   - Performs MSA using MAFFT with FFT-NS-2 mode for speed and auto-configured threading.
   - Dynamically determines the optimal number of CPU threads based on system resources.
   - Usage:
   
   ```bash
   python scripts/run_alignment.py
   ```

---

## Running the Pipeline

1. **Prepare the Raw Data:**

   - Place your viral genome sequences (in FASTA format) in the `data/raw/`directory.

2. **Run Preprocessing:**

   - Execute the preprocessing script:
   ```bash
   python scripts/preprocess_data.py
   ```
   - This generates cleaned_sequences.fasta in the `data/processed/` directory.

3. **Run Alignment and Conserved Region Analysis:**

   - Execute the alignment script:
   ```bash
   python scripts/run_alignment.py
   ```
   - This generates:
      - `aligned_sequences.fasta` (aligned sequences).
      - `conserved_regions.npy` (conserved positions).

4. **Explore the Notebooks:**
   
   - Use the notebooks for additional visualization and experimentation. 

---

## Next Steps

The Cas13Hunter pipeline is a work in progress, with several planned enhancements to improve its functionality, scalability, and accuracy. Here are the key areas of development:

1. **Off-Target Analysis and Validation of Guide RNAs**:
   - Incorporate off-target analysis to ensure designed guide RNAs are specific and minimize unintended targeting.
   - Validate guide RNAs against reference genomes and other viral sequences to identify and filter potential off-target effects.

2. **Additional Preprocessing Steps**:
   - Expand the current preprocessing workflow beyond filtering, deduplication, and GC content calculations.
   - Implement advanced quality control and sequence filtering methods to handle complex datasets and ensure high-quality input for downstream analyses.

3. **Enhanced MSA Solutions**:
   - Introduce state-of-the-art alternatives to MAFFT for multiple sequence alignment, such as:
     - **Clustal Omega** and **MUSCLE** for improved accuracy.
     - **Kalign3**, a GPU-accelerated tool for faster alignment.
   - Evaluate these alternatives to improve both alignment speed and accuracy.

4. **Scalability**:
   - Optimise the pipeline to handle large-scale datasets efficiently.
   - Scale the pipeline to process thousands or even millions of sequences through:
     - Code optimisation.
     - Leveraging HPC cluster for testing and performance enhancement.

These enhancements will ensure that Cas13Hunter remains a robust and scalable tool for CRISPR-Cas13 target identification and RNA therapeutic design.


---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.


---

## Contact

For questions or collaboration inquiries, please reach out to `hsellak@outlook.com`.

> **Note:** Contributions are welcome! Fork the repository, make your changes, and submit a pull request.