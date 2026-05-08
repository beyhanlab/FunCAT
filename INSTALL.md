# FunCAT Installation

**Fungal Chromosome Assembly Tool**

---

## One-command install (recommended)

```bash
git clone https://github.com/beyhanlab/FunCAT.git
cd FunCAT
conda env create -f environment.yml
conda activate funcat
funcat
```

That's it. Everything installs in one step — Flye, Medaka, Polypolish, Racon,
BWA, samtools, bcftools, bgzip, tabix, purge_dups, MUMmer, and FunCAT itself.

---

## What gets installed

### Core tools (required)
| Tool | Purpose |
|------|---------|
| flye | De novo assembly |
| minimap2 | Read mapping |
| seqkit | Read statistics |
| filtlong | Read length filtering |
| samtools ≥1.17 | BAM handling and coverage |
| bcftools | Variant calling (required by Medaka) |
| htslib (bgzip + tabix) | VCF compression (required by Medaka) |
| medaka | Nanopore consensus polishing |

### Optional tools (unlock additional modules)
| Tool | Module unlocked |
|------|----------------|
| racon | PacBio polishing |
| bwa | Illumina hybrid polishing |
| polypolish + polypolish_insert_filter | Illumina hybrid polishing |
| purge_dups (+ split_fa, pbcstat, calcuts, get_seqs) | Diploid haplotig removal |
| mummer4 (nucmer) | SNP analysis and dotplots |
| pilon | Alternative Illumina polisher |

> **Note:** Optional tools are still installed by `environment.yml`.
> FunCAT will warn you gracefully if any are missing rather than crashing.

---

## Manual install

If you prefer to manage your own environment:

```bash
# 1. Create environment
conda create -n funcat python=3.10
conda activate funcat

# 2. Core required tools
conda install -c bioconda -c conda-forge \
  flye minimap2 seqkit filtlong \
  "samtools>=1.17" bcftools htslib

# 3. Medaka (pip only — no stable conda package)
#    bcftools + htslib must be installed BEFORE medaka
pip install medaka pyabpoa

# 4. Optional tools
conda install -c bioconda racon purge_dups mummer4 bwa polypolish

# 5. Install FunCAT
pip install -e .
```

---

## Common install issues

**Medaka crashes with "bcftools not found" or "tabix not found"**
```bash
conda install -c bioconda bcftools htslib
```

**Old samtools (0.1.x) conflict**
If you see `samtools: invalid option -- 'a'` or similar, you have an ancient
samtools on your PATH. Fix with:
```bash
conda install -c bioconda "samtools>=1.17"
which samtools   # should now point to your conda env
```

**polypolish_insert_filter not found**
This binary is included in the polypolish conda package. Reinstall with:
```bash
conda install -c bioconda polypolish
```

---

## Verify installation

```bash
funcat --help
```

## Requirements

- Python 3.9+
- conda (strongly recommended — most tools are not pip-installable)
- macOS or Linux (Windows not supported)
