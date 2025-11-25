#!/bin/bash
set -euo pipefail

#######
###! Genome-wide association workflow; ~250 samples recommended
#######

# Variable definitions for easier reuse

INPUT_DIR="example"
SCRIPT_DIR="12-GWAS/script"
MERGE_FASTA_DIR="merge_fasta"
PYTHON3="python3"
VERYFASTTREE="VeryFastTree"

# 1. Extract sequences for the target population
python3 "$SCRIPT_DIR/0_cancat.py" \
    "$INPUT_DIR/list.txt" \
    "$MERGE_FASTA_DIR" \
    "$INPUT_DIR/alignmentSub.aln"

# 2. Extract SNP positions
snp-sites -c -o "$INPUT_DIR/alignmentSub.temp" "$INPUT_DIR/alignmentSub.aln"

# 3. Build the phylogenetic tree
nohup "$VERYFASTTREE" -nt -threads 16 \
    "$INPUT_DIR/alignmentSub.temp" \
    >"$INPUT_DIR/treeFull.nwk" \
    2>"$INPUT_DIR/treeFull.log" &

# 4. Generate genotype file (VCF)
snp-sites -c -v -o "$INPUT_DIR/tmp1.vcf" "$INPUT_DIR/alignmentSub.aln"

# 5. Process the VCF file and extract required columns
sed "s/#CHROM/CHROM/" "$INPUT_DIR/tmp1.vcf" \
    | grep -v "#" | cut -f 2,4-5,10- | sed "1s/POS/ps/" \
    >"$INPUT_DIR/tmp2.txt"

# 6. Create the bugwas input file (biallelic filter enabled by default)
${PYTHON3} 0_prepare.py \
    "$INPUT_DIR/tmp2.txt" \
    --output_file "$INPUT_DIR/geno_biallelic_SNP.txt"


# 7. Run GWAS analysis and post-processing (requires phenotype and metadata files)
Rscript s5_analysis_GWAS.r

# Gene annotation example
${PYTHON3} 12-GWAS/script/3_annotation.py \
    --gwas_file "bugwas_biallelic_lmmout_allSNPs.txt" \
    --gff_file "NC_000915.gff" \
    --output_dir "output" \
    --logp_threshold 5 \
    --dist_threshold 1

# 8. Remove intermediate files (comment out to keep them)
rm -f "$INPUT_DIR/tmp1.vcf" "$INPUT_DIR/alignmentSub.aln" "$INPUT_DIR/alignmentSub.temp"

