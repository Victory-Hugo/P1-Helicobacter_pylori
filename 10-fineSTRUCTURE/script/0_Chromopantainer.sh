#!/bin/bash
# Input files are BED/BIM/FAM converted from a VCF
# Use PLINK to convert VCF to BED format
#!======Required software before running=======
#todo 1. PLINK
#todo 2. Shapeit
#todo 3. Perl
#todo 4. ChromoPainterv2
#todo 5. fineSTRUCTURE
#todo 6. R (for visualization)
#!=================================
# Resolve project directories relative to this script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
INPUT_DIR="${BASE_DIR}/data" # Provide vcf.gz and its index here
OUTPUT_DIR="${BASE_DIR}/OUTPUT/CDS" # Should be empty before running
mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"

#todo Configuration files located in ./conf/
GENETIC_MAP_FINAL="${BASE_DIR}/conf/genetic_map_Finalversion_HP.txt" # Derived from genetic_map_HP.txt
GENETIC_MAP="${BASE_DIR}/conf/genetic_map_HP.txt" # HP genetic map computed previously
D_R_FILE="${BASE_DIR}/conf/popDonRec.txt"

# Define input VCF and output prefixes
VCF_FILE="${INPUT_DIR}/merged.vcf.gz"
PLINK_PREFIX="${OUTPUT_DIR}/CDS"
PLINK_HP_PREFIX="${PLINK_PREFIX}_HP"


# # Step 1: Filter VCF with PLINK and generate BED format
plink --vcf "$VCF_FILE" \
      --make-bed \
      --double-id \
      --allow-extra-chr \
      --out "$PLINK_PREFIX" \
      --geno 0.1 \
      --mind 0.1 \
      --indep-pairwise 50 10 0.1 

# Backup PLINK outputs
FILES=("bed" "bim" "fam")
for EXT in "${FILES[@]}"; do
    FILE="${PLINK_PREFIX}.${EXT}"
    cp "$FILE" "${FILE}.bak"
done

# Step 2: Replace chromosome identifiers in the BIM file
awk 'BEGIN{OFS="\t"} { $1 = ($1 == "NC_000915.1" ? "1" : $1); print }' \
     "${PLINK_PREFIX}.bim" > \
     "${OUTPUT_DIR}/temp.bim" && mv "${OUTPUT_DIR}/temp.bim" "${PLINK_PREFIX}.bim"

echo "Replacement finished. The following file was updated:"
echo "${PLINK_PREFIX}.bim"
echo "Original file backed up at: ${PLINK_PREFIX}.bim.bak"

# Step 3: Use PLINK to extract chromosome 1 and create a new BED set
plink --bfile "$PLINK_PREFIX" \
      --chr 1 \
      --make-bed \
      --alleleACGT \
      --out "$PLINK_HP_PREFIX"

# Step 4: Phase with Shapeit
shapeit --input-bed "${PLINK_HP_PREFIX}.bed" "${PLINK_HP_PREFIX}.bim" "${PLINK_HP_PREFIX}.fam" \
             --input-map "$GENETIC_MAP" \
             --output-max "${PLINK_PREFIX}_HP.phased.haps" "${PLINK_HP_PREFIX}.phased.sample" \
             --output-log "${PLINK_HP_PREFIX}.log" \
             --force \
             --burn 10 \
             --prune 10 \
             --main 30 \
             --thread 16

# # Step 5: Generate ID file for ChromoPainter
# #! Note: Manually replace the second column with group names and the third column with 1 after this step.
awk 'BEGIN{FS=" "; OFS="\t"} {print $2, $1, $6}' "${PLINK_HP_PREFIX}.fam" | tr -s '\t ' ' ' > \
     "${PLINK_HP_PREFIX}.ids"


# Step 6: Convert Shapeit output to fineSTRUCTURE format using Perl helpers
echo "Converting via impute2; this step may take a long time."
impute2chromopainter.pl -J \
     "${PLINK_HP_PREFIX}.phased.haps" \
     "${PLINK_HP_PREFIX}.phase"
dos2unix "$GENETIC_MAP_FINAL"
echo 'Ensure genetic_map_Finalversion_HP.txt is in UNIX format.'



convertrecfile.pl -M \
     hap "${PLINK_HP_PREFIX}.phase" \
     "$GENETIC_MAP_FINAL" \
     "${PLINK_HP_PREFIX}.recombfile"

echo "All steps completed."
# Remove backup .bak files
rm "${PLINK_PREFIX}.bim.bak"
rm "${PLINK_PREFIX}.bed.bak"
rm "${PLINK_PREFIX}.fam.bak"
echo "Backup files removed."

dos2unix "${PLINK_HP_PREFIX}.ids"
echo "Ensured ${PLINK_HP_PREFIX}.ids is in UNIX format."
dos2unix "${D_R_FILE}"
echo "Ensured ${D_R_FILE} is in UNIX format."
