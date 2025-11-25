#!/usr/bin/env bash
set -euo pipefail

# Input VCF file
vcf="WGS.vcf.gz"
# Output prefix (customize as needed)
out_prefix="subset_WGS"

# (1) Use PLINK to generate binary files
plink --vcf "$vcf" \
      --double-id \
      --make-bed \
      --allow-extra-chr \
      --out "$out_prefix"

# Back up the original .bim
cp "${out_prefix}.bim" "${out_prefix}.bim.bak"

# (2) Fix chromosome name: replace NC_000915.1 with 1
awk 'BEGIN{OFS="\t"} 
     $1=="NC_000915.1"{$1="1"} 
     {print}' \
    "${out_prefix}.bim" > "${out_prefix}.tmp.bim" \
&& mv "${out_prefix}.tmp.bim" "${out_prefix}.bim"

echo "✅ ${out_prefix}.bim has been updated; original file saved as ${out_prefix}.bim.bak"

# (3) Optional: filter loci/samples with high missingness
# Remove --alleleACGT for PLINK1.9; keep it for PLINK2.0
plink --bfile "$out_prefix" \
      --geno 0.10 \
      --mind 0.10 \
      --indep-pairwise 50 10 0.1 \
      --max-maf 0.01 \
      --make-bed \
      --alleleACGT \
      --out "${out_prefix}_filtered"

echo "✅ Filtering complete; results stored in ${out_prefix}_filtered.{bed,bim,fam}"
