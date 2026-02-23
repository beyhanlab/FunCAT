# FungalFlye

Long-read fungal genome assembly pipeline using Flye.

## Install

### Option 1 — Conda (recommended)

conda env create -f environment.yml  
conda activate fungalflye  

pip install git+https://github.com/beyhanlab/FungalFlye.git

## Run

fungalflye

## Features

- Read length filtering
- Coverage normalization
- Flye assembly
- Racon polishing
- Redundant contig pruning
- Telomere completeness detection
- Assembly QC report

## Requirements

Nanopore long reads (FASTQ)
