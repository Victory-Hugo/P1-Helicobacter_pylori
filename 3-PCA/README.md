# PCA Analysis for Helicobacter pylori

This directory contains scripts to perform principal component analysis (PCA) on *Helicobacter pylori* variants and to visualize the PCA results.

## Directory structure

- `script/0-Plink_PCA.sh`  
  Run PLINK (v1.9) to compute PCA from a VCF file.
- `script/1-PCA_vis.R`  
  Visualize PCA results as a scatter plot PDF.
- `script/conf/` (optional)  
  Place configuration files or helper scripts here if needed.

## 1. PCA with PLINK

Script: `script/0-Plink_PCA.sh`

### Input

- `../cleaned.vcf.gz`  
  Compressed VCF file with filtered variants.

### Output

- `../filtered/`  
  PLINK binary files (`filtered.bed`, `filtered.bim`, `filtered.fam`, LD-pruned files, etc.).
- `../filtered_PCA/hpglobal_LD_PCA.eigenval`  
  Eigenvalues from PCA.
- `../filtered_PCA/hpglobal_LD_PCA.eigenvec`  
  Eigenvectors (principal components) per sample.

### Usage

From the `3-PCA` directory:

```bash
bash script/0-Plink_PCA.sh
```

Requirements:

- PLINK v1.9 (or compatible) available in `PATH`.
- Input VCF at `../cleaned.vcf.gz`.

## 2. PCA visualization in R

Script: `script/1-PCA_vis.R`

This script reads a prepared PCA table and a color mapping file, and produces a PDF scatter plot.

### Required input files

All files are read from the current working directory when running the script.

- `PCA.csv`  
  Must contain at least the following columns:
  - `ID`
  - `PC1`
  - `PC2`
  - `Class_big`
  - `Class_small`

- `color.csv`  
  Must contain:
  - `Class_big`
  - `color` (any valid R color string or hex code, e.g. `#1f77b4`)

### Output

- `PCA.pdf`  
  A PCA scatter plot, colored by `Class_big` and shaped by `Class_small`.  
  The script writes console messages describing the mapping and number of points.

### Usage

From the directory that contains `PCA.csv` and `color.csv`:

```bash
Rscript script/1-PCA_vis.R
```

Or in an interactive R session:

```r
source("script/1-PCA_vis.R")
```

Required R packages:

- `tidyverse`
- `RColorBrewer`

## Notes

- Paths in the scripts are relative so they can be moved as a group without exposing machine-specific absolute paths.
- If you change where the input VCF or PCA tables are stored, update the relative paths in the scripts accordingly.

