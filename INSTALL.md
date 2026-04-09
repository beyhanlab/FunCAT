# FunCAT Installation

**Fungal Chromosome Assembly Tool**

## One-command install (recommended)

```bash
git clone https://github.com/beyhanlab/FunCAT.git
cd FunCAT
conda env create -f environment.yml
conda activate funcat
funcat
```

That's it. Everything — Flye, Medaka, Polypolish, Racon, BWA, samtools, purge_dups, and FunCAT itself — installs in one step.

## Manual install

If you prefer to set up your own environment:

```bash
# Create environment
conda create -n funcat python=3.10
conda activate funcat

# Required tools
conda install -c bioconda flye minimap2 seqkit filtlong samtools

# ONT polisher
pip install medaka

# Optional tools (unlock additional modules)
conda install -c bioconda racon purge_dups mummer bwa
pip install polypolish

# Install FunCAT
pip install -e .
```

## Verify installation

```bash
funcat --help
```

## Requirements

- Python 3.9+
- conda (recommended) or pip
- macOS or Linux
