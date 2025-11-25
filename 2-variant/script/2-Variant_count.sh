#!/bin/bash

# Set working directories based on script location
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" &> /dev/null && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
DATA_DIR="$PROJECT_DIR/data"
OUTPUT_DIR="$PROJECT_DIR/output"
TEMP_DIR="$OUTPUT_DIR/temp"

# Ensure output directories exist
mkdir -p "$OUTPUT_DIR"
mkdir -p "$TEMP_DIR"

# Set log file
LOG_FILE="$OUTPUT_DIR/variant_stats.log"
echo "Starting variant statistics analysis: $(date)" > "$LOG_FILE"

# Define frequency classification criteria
# Common variants: MAF > 0.05 (5%)
# Low-frequency variants: 0.01 < MAF ≤ 0.05 (1–5%)
# Rare variants: MAF ≤ 0.01 (≤1%)
# Singleton: variant observed in only one sample
# Doubleton: variant observed in exactly two samples

# Create summary result file
SUMMARY_FILE="$OUTPUT_DIR/variant_summary.tsv"
echo -e "Region\tTotal_Variants\tCommon\tLow_Frequency\tRare\tSingleton\tDoubleton\tPrivate" > "$SUMMARY_FILE"

# Create directory for variant IDs used in Venn diagrams
VENN_DIR="$OUTPUT_DIR/venn_data"
mkdir -p "$VENN_DIR"

# Process each continental VCF file
for VCF_FILE in "$DATA_DIR"/*.vcf.gz; do
    # Use file name as region name
    REGION=$(basename "$VCF_FILE" .vcf.gz)
    echo "Processing $REGION... $(date)" >> "$LOG_FILE"
    
    # 1. Count total variants
    TOTAL_VARIANTS=$(bcftools view -H "$VCF_FILE" | wc -l)
    echo "  Total variants: $TOTAL_VARIANTS" >> "$LOG_FILE"
    
    # 2. Generate allele frequency information
    echo "  Calculating allele frequencies..." >> "$LOG_FILE"
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AC\t%AN\n' "$VCF_FILE" > "$TEMP_DIR/${REGION}_freq.txt"
    
    # 3. Classify variants based on frequency
    # Add MAF column and classify
    echo "  Classifying variants by frequency..." >> "$LOG_FILE"
    
    # Use AWK to process frequency data
    awk -v region="$REGION" '
    BEGIN {
        common = 0;
        low_freq = 0;
        rare = 0;
        singleton = 0;
        doubleton = 0;
        
        
        common_file = "'$VENN_DIR'/" region "_common.txt";
        low_file = "'$VENN_DIR'/" region "_low.txt";
        rare_file = "'$VENN_DIR'/" region "_rare.txt";
        singleton_file = "'$VENN_DIR'/" region "_singleton.txt";
        doubleton_file = "'$VENN_DIR'/" region "_doubleton.txt";
        all_variants_file = "'$VENN_DIR'/" region "_all.txt";
    }
    {
        # Calculate MAF (Minor Allele Frequency)
        ac = $5;
        an = $6;
        if (an == 0) next;  # skip sites without genotype information
        
        # Calculate allele frequency
        af = ac / an;
        maf = (af > 0.5) ? (1 - af) : af;  # take the smaller allele frequency
        
        # Variant ID (used for Venn diagram)
        var_id = $1 "_" $2 "_" $3 "_" $4;
        print var_id >> all_variants_file;
        
        # Classify by frequency
        if (ac == 1) {
            singleton++;
            print var_id >> singleton_file;
        } else if (ac == 2) {
            doubleton++;
            print var_id >> doubleton_file;
        }
        
        if (maf > 0.05) {
            common++;
            print var_id >> common_file;
        } else if (maf > 0.01 && maf <= 0.05) {
            low_freq++;
            print var_id >> low_file;
        } else {
            rare++;
            print var_id >> rare_file;
        }
    }
    END {
        print "Common variants (MAF > 0.05): " common > "'$LOG_FILE'";
        print "Low-frequency variants (0.01 < MAF ≤ 0.05): " low_freq >> "'$LOG_FILE'";
        print "Rare variants (MAF ≤ 0.01): " rare >> "'$LOG_FILE'";
        print "Singleton (AC=1): " singleton >> "'$LOG_FILE'";
        print "Doubleton (AC=2): " doubleton >> "'$LOG_FILE'";
        
        # Save classification summary data
        print "'$REGION'\t" (common + low_freq + rare) "\t" common "\t" low_freq "\t" rare "\t" singleton "\t" doubleton "\t0" > "'$TEMP_DIR/${REGION}_summary.txt'";
    }' "$TEMP_DIR/${REGION}_freq.txt"
    
    # Append results to summary file
    cat "$TEMP_DIR/${REGION}_summary.txt" >> "$SUMMARY_FILE"
    
    echo "Finished processing $REGION" >> "$LOG_FILE"
done

# 4. Calculate private variants - find variants that appear only in a single region
echo "Calculating private variants..." >> "$LOG_FILE"

ALL_VARIANTS_DIR="$TEMP_DIR/all_variants"
mkdir -p "$ALL_VARIANTS_DIR"

for REGION_FILE in "$VENN_DIR"/*_all.txt; do
    REGION=$(basename "$REGION_FILE" _all.txt)
    cp "$REGION_FILE" "$ALL_VARIANTS_DIR/${REGION}.txt"
done

# Calculate private variants for each region
for REGION_FILE in "$ALL_VARIANTS_DIR"/*.txt; do
    REGION=$(basename "$REGION_FILE" .txt)
    
    # Create a file containing variants from all other regions
    find "$ALL_VARIANTS_DIR" -type f ! -name "${REGION}.txt" -exec cat {} \; | sort | uniq > "$TEMP_DIR/other_regions.txt"
    
    # Identify variants present only in the current region (private variants)
    comm -23 <(sort "$REGION_FILE" | uniq) <(sort "$TEMP_DIR/other_regions.txt") > "$VENN_DIR/${REGION}_private.txt"
    
    # Count private variants
    PRIVATE_COUNT=$(wc -l < "$VENN_DIR/${REGION}_private.txt")
    
    # Update private variant counts in the summary file
    awk -v region="$REGION" -v private="$PRIVATE_COUNT" '
        $1 == region {$7 = private; print $0; next}
        {print}
    ' "$SUMMARY_FILE" > "$TEMP_DIR/summary_temp.tsv" && mv "$TEMP_DIR/summary_temp.tsv" "$SUMMARY_FILE"
    
    echo "$REGION private variants: $PRIVATE_COUNT" >> "$LOG_FILE"
done

# Generate files suitable for Venn diagram analysis
echo "Generating Venn diagram data files..." >> "$LOG_FILE"

# Reformat summary data for easier visualization
cat "$SUMMARY_FILE" | column -t > "$OUTPUT_DIR/variant_summary_formatted.tsv"
