#!/bin/bash

# Get the absolute path of the script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" &> /dev/null && pwd)"
# Get the project directory from the script directory
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
# Get the root directory (parent of the project directory)
ROOT_DIR="$(dirname "$PROJECT_DIR")"
# Set data output directory
DATA_DIR="$PROJECT_DIR/data"
# Set configuration directory (sample lists)
CONF_DIR="$PROJECT_DIR/conf"
# Set path to merged VCF file (relative to project root)
MERGED_VCF="$ROOT_DIR/1-nucmer/output/merge/merged_biallelic.7544.vcf.gz"

# Ensure output directory exists
mkdir -p "$DATA_DIR"

# Use GNU parallel for parallel processing
find "$CONF_DIR" -name "*.txt" | \
parallel -j 8 'bcftools view -S {} -Oz -o '"$DATA_DIR"'/$(basename {} .txt).vcf.gz '"$MERGED_VCF"

echo "All VCF files have been generated!"
