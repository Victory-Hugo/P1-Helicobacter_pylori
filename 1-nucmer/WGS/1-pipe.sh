gatk VariantFiltration \
    -R "${REF}" \
    -V "${RAW_VCF}" \
    -O "${FILTERED_FLAGGED_VCF}" \
    --filter-name "LOW_QUAL" \
    --filter-expression "QUAL < 20.0 || QD < 1.0 || MQ < 30.0 || FS > 100.0 || SOR > 10.0"

vcftools \
    --vcf "${FILTERED_FLAGGED_VCF}" \
    --remove-filtered-all \
    --recode \
    --recode-INFO-all \
    --out "${OUT_DIR}/step1_pass"

mv "${OUT_DIR}/step1_pass.recode.vcf" "${FILTERED_STEP1_VCF}"

FILTERED_STEP2_VCF="${OUT_DIR}/step2_highconf.recode.vcf"

vcftools \
    --vcf "${FILTERED_STEP1_VCF}" \
    --max-missing 0.95 \
    --maf 0.05 \
    --minDP 3 \
    --min-meanDP 20 \
    --minQ 30 \
    --recode \
    --recode-INFO-all \
    --out "${OUT_DIR}/step2_highconf"

: > "${LOW_QUAL_SAMPLES}"

REF_LEN=$(awk '{if($0 !~ /^>/) len += length($0)} END{print len}' "${REF}")

while read -r SAMPLE; do
    [ -z "${SAMPLE}" ] && continue

    ASM="${ASM_DIR}/${SAMPLE}.fa"
    BAM="${BAM_DIR}/${SAMPLE}.bam"

    if [ ! -f "${ASM}" ]; then
        :
    else
        N_CONTIGS=$(grep -c "^>" "${ASM}")
        if [ "${N_CONTIGS}" -gt 500 ]; then
            echo "${SAMPLE}" >> "${LOW_QUAL_SAMPLES}"
            continue
        fi
    fi

    if [ ! -f "${BAM}" ]; then
        continue
    fi

    if [ ! -f "${BAM}.bai" ]; then
        samtools index "${BAM}"
    fi

    COVERED_BASES=$(samtools depth -a "${BAM}" | awk '{if($3>0) cov++} END{print (cov>0?cov:0)}')
    if [ "${REF_LEN}" -eq 0 ]; then
        exit 1
    fi

    COV_PCT=$(awk -v cov="${COVERED_BASES}" -v total="${REF_LEN}" 'BEGIN{printf "%.2f", cov/total*100}')

    BELOW=$(awk -v x="${COV_PCT}" 'BEGIN{print (x<70)?1:0}')
    if [ "${BELOW}" -eq 1 ]; then
        echo "${SAMPLE}" >> "${LOW_QUAL_SAMPLES}"
    fi

done < "${SAMPLES_LIST}"

if [ -s "${LOW_QUAL_SAMPLES}" ]; then
    vcftools \
        --vcf "${FILTERED_STEP2_VCF}" \
        --remove "${LOW_QUAL_SAMPLES}" \
        --recode \
        --recode-INFO-all \
        --out "${OUT_DIR}/final_filtered"
else
    cp "${FILTERED_STEP2_VCF}" "${OUT_DIR}/final_filtered.recode.vcf"
fi
