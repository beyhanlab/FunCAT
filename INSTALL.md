# FungalFlye — Installation

## Fresh install

```bash
git clone https://github.com/beyhanlab/FungalFlye.git
cd FungalFlye
conda create -n fungalflye python=3.10 -y
conda activate fungalflye
conda install -c bioconda flye minimap2 seqkit filtlong samtools racon -y
pip install medaka
pip install -e .
fungalflye --help
```

## Updating (if already installed)

```bash
cd /path/to/FungalFlye
conda activate fungalflye
git pull origin main
```

That's it. No reinstall needed.

## Notes

- Always use `pip install medaka` — not conda
- Always activate the `fungalflye` environment before running the tool
- Run `fungalflye` from any directory once installed
