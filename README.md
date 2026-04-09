# FunCAT 🧬🍄

**Fungal Chromosome Assembly Tool**

Long-read fungal genome assembly from raw Nanopore or PacBio HiFi reads to a polished, chromosome-scale assembly — in a single command.

Developed by Jacob Durazo · Beyhan Lab · J. Craig Venter Institute

---

## Install

```bash
git clone https://github.com/beyhanlab/FunCAT.git
cd FunCAT
conda env create -f environment.yml
conda activate funcat
```

> **Don't have conda?** Install [Miniconda](https://docs.conda.io/en/latest/miniconda.html) first.

## Run

```bash
funcat
```

That's it. The interactive wizard guides you through everything.

---

## What it does

FunCAT takes raw Nanopore or PacBio HiFi reads and produces chromosome-scale fungal genome assemblies through a fully automated pipeline:

- **Adaptive Flye assembly** — auto-tunes parameters from your read N50, coverage, and GC content
- **Iterative Medaka polishing** — polishes up to 3 rounds, stops when the assembly converges
- **Repeat-aware scaffolding** — uses long reads to bridge contig gaps at repeat boundaries
- **Telomere-guided scaffolding** — attaches telomeric fragments to uncapped chromosome ends
- **Illumina hybrid polishing** — optional short-read polishing to reduce internal stop codons
- **Contig confidence scoring** — flags collapsed repeats and contamination before you publish
- **Auto-generated clean assembly** — removes flagged contigs, produces publication-ready FASTA
- **Self-contained HTML report** — share a single file with collaborators, no internet required

---

## Benchmark

Tested on *Histoplasma capsulatum* Nanopore HQ reads (18x coverage):

| | Vanilla Flye | FunCAT |
|---|---|---|
| Contigs | 219 | 36 |
| N50 | 1.6 Mb | 3.34 Mb |
| Largest contig | 2.65 Mb | 6.06 Mb |
| Runtime | ~45 min | ~65 min |

**83% reduction in contig count. N50 more than doubled. Same raw data.**

---

## Read types

| Flag | Reads | Polisher |
|---|---|---|
| `nano-hq` | Nanopore R10.4+, Q20 (default) | Medaka |
| `nano-raw` | Nanopore R9.4, standard accuracy | Medaka |
| `pacbio-hifi` | PacBio HiFi / CCS | Racon |

---

## Output files

```
outdir/
├── final.fasta                  # Full assembly (all contigs)
├── final_clean.fasta            # Clean assembly (FLAG contigs removed)
├── nuclear.fasta                # Nuclear contigs only
├── mitochondrial.fasta          # Mitochondrial contigs (if detected)
├── final_funcat_report.html     # Self-contained HTML report
├── confidence/
│   └── contig_confidence.tsv   # Per-contig coverage scores
└── fungalflye_qc/
    ├── stats.txt
    ├── telomeres.tsv
    └── length_histogram.png
```

---

## Citation

If you use FunCAT in published work, please cite:

> Durazo J, et al. FunCAT: a purpose-built long-read assembly toolkit for fungal genomes. (2025) [Manuscript in preparation]

Also cite the underlying tools:
- **Flye:** Kolmogorov et al. (2019) *Nature Biotechnology*
- **Medaka:** Oxford Nanopore Technologies
- **minimap2:** Li (2018) *Bioinformatics*
- **purge_dups:** Guan et al. (2020) *Genome Biology*

---

## License

MIT License — see [LICENSE](LICENSE) for details.

*Developed by the Beyhan Lab, J. Craig Venter Institute*
