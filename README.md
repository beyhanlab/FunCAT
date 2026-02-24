# FungalFlye

Long-read fungal genome assembly pipeline using Flye.


⭐ Option 1 — Conda (Recommended)

This installs all dependencies + FungalFlye in one environment.

conda env create -f environment.yml
conda activate fungalflye

pip install git+https://github.com/beyhanlab/FungalFlye.git


▶️ Run FungalFlye

Launch the interactive wizard:

fungalflye

That’s it. The pipeline will guide you through the rest.


🧬 What FungalFlye Does

FungalFlye is an end-to-end long-read fungal genome assembly pipeline designed for chromosome-level assemblies with minimal effort.

Features include:

✨ Intelligent read length filtering
✨ Optional coverage normalization / downsampling
✨ Flye assembly with fungal-optimized parameters
✨ Racon polishing
✨ Automatic redundant contig pruning (>95% containment)
✨ Telomere motif discovery and chromosome completeness detection
✨ Publication-ready assembly QC metrics
✨ Interactive step-by-step wizard interface

📊 Outputs

FungalFlye produces:

Final polished assembly FASTA

Contig size tables

Assembly statistics (N50, L50, genome size)

Telomere completeness report

QC plots and summaries

All results are organized automatically in the output directory.


📦 Requirements

Input:

Oxford Nanopore long reads (FASTQ or FASTQ.GZ)

Recommended:

≥40× genome coverage

High-molecular-weight DNA
