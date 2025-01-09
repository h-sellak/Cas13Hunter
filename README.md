# Cas13Hunter

A bioinformatics pipeline for identifying conserved RNA regions in viral genomes and designing CRISPR-Cas13 guide RNAs for therapeutic development. This tool focuses on enabling rapid target identification for RNA respiratory viruses like SARS-CoV-2, leveraging multiple sequence alignment, functional annotation, and target prediction techniques.

## Getting Started

These instructions will help you set up the project environment and get started with Cas13Hunter.

---

## Prerequisites

To run this project, ensure you have the following installed:

1. [Anaconda/Miniconda](https://docs.conda.io/en/latest/miniconda.html)
2. [Git](https://git-scm.com/)

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

To maintain consistency, here is the recommended directory structure for the project:

```
Cas13Hunter/
├── data/                 # Raw and processed viral genome data
├── src/                  # Python scripts for pipeline modules
├── notebooks/            # Jupyter notebooks for experimentation
├── tests/                # Unit tests for the pipeline
├── environment.yml       # Conda environment configuration file
└── README.md             # Project documentation
```

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
