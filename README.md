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
├── scripts/               # Python scripts for pipeline modules
│   ├── preprocess_data.py             # Script for sequence preprocessing
│   ├── run_alignment.py               # Script for running MAFFT alignment
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

1. **Data Preparation**:
   - Collect viral genome sequences in FASTA/FASTQ format and save them in the `data/raw/` directory.

2. **Run the Pipeline**:
   - Coming soon: Detailed instructions on running the sequence alignment and analysis pipeline.

3. **Contribute**:
   - Contributions are welcome! Fork the repository, make your changes, and submit a pull request.

---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## Contact

For questions or collaboration inquiries, please reach out to `hsellak@outlook.com`.
