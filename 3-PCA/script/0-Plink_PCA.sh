#!/bin/bash

# Script name: simple_pca.sh
# Function: Perform principal component analysis (PCA) on a VCF file using PLINK (v1.9)


# Error handling
set -e

# Define input and output paths (relative paths)
INPUT_VCF="../cleaned.vcf.gz"
OUTPUT_DIR="../filtered"
mkdir -p "$OUTPUT_DIR"
PLINK_PREFIX="${OUTPUT_DIR}/filtered"
PCA_OUTPUT_DIR="../filtered_PCA"
mkdir -p "$PCA_OUTPUT_DIR"

# Convert VCF to PLINK binary format
echo "Converting VCF to PLINK binary format"
plink --vcf "$INPUT_VCF" --make-bed \
      --double-id  \
      --allow-extra-chr \
      --out "$PLINK_PREFIX"

# Perform linkage disequilibrium (LD) pruning
echo "Performing linkage disequilibrium (LD) pruning"
plink --bfile "$PLINK_PREFIX" \
      --indep-pairwise 50 10 0.1 \
	--allow-extra-chr \
      --out "${PLINK_PREFIX}_ld_pruned"

# Run principal component analysis (PCA)
echo "Running principal component analysis (PCA)"
plink --bfile "$PLINK_PREFIX" \
      --extract "${PLINK_PREFIX}_ld_pruned.prune.in" \
      --pca 20 \
	  --allow-extra-chr \
      --out "${PCA_OUTPUT_DIR}/hpglobal_LD_PCA"

echo "PCA analysis completed! All output files are saved in ${PCA_OUTPUT_DIR}/"
