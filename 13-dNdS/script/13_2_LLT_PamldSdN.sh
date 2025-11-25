#!/bin/bash
##############
## 1. dNdS ##
##############

# Exit on errors and report unset variables
set -euo pipefail

LIST_DIR="list"
OUTPUT_ROOT="dSdN/OUTPUT"
YN00_EXECUTABLE="paml-master/bin/yn00"
export LIST_DIR OUTPUT_ROOT YN00_EXECUTABLE

# Define worker that processes each ID
process_id() {
    local i="$1"
    # Build file paths
    COMPARISONS_FILE="${LIST_DIR}/comparisons_CDS_${i}_core.txt"
    OUTPUT_BASE_DIR="${OUTPUT_ROOT}/CDS_${i}"
    YN00_EXEC="${YN00_EXECUTABLE}"
    CONTROL_FILE="yn00.ctl"

    # Print an error and exit
    error_exit() {
        echo "Error: $1" >&2
        exit 1
    }

    # Ensure comparisons file exists
    if [[ ! -f "$COMPARISONS_FILE" ]]; then
        error_exit "Comparisons file $COMPARISONS_FILE not found."
    fi

    # Ensure yn00 is available
    if [[ ! -x "$YN00_EXEC" ]]; then
        error_exit "yn00 executable $YN00_EXEC does not exist or is not executable."
    fi

    # Ensure control file exists
    if [[ ! -f "$CONTROL_FILE" ]]; then
        error_exit "Control file $CONTROL_FILE not found."
    fi

    # Loop through each row in comparisons.txt
    while IFS=$'\t' read -r ALN1 ALN2 || [[ -n "${ALN1:-}" && -n "${ALN2:-}" ]]; do
        # Skip empty or malformed rows
        if [[ -z "$ALN1" || -z "$ALN2" ]]; then
            echo "Skipping empty or malformed row."
            continue
        fi

        echo "Processing files: $ALN1 and $ALN2"

        # Check the two .aln files
        if [[ ! -f "$ALN1" ]]; then
            echo "Warning: File $ALN1 not found, skipping this pair."
            continue
        fi

        if [[ ! -f "$ALN2" ]]; then
            echo "Warning: File $ALN2 not found, skipping this pair."
            continue
        fi

        # Derive base names (without path or extension)
        BASENAME1=$(basename "$ALN1" .aln)
        BASENAME2=$(basename "$ALN2" .aln)

        # Build new output directory name
        NEW_OUTPUT_DIR="${OUTPUT_BASE_DIR}/${BASENAME1}_${BASENAME2}"

        # Ensure the output directory exists
        mkdir -p "$NEW_OUTPUT_DIR"
        echo "Created output directory: $NEW_OUTPUT_DIR"

        # Merge the two .aln files into NEW_OUTPUT_DIR/temp.aln
        MERGE_ALN="${NEW_OUTPUT_DIR}/temp.aln"
        {
            cat "$ALN1"
            echo ""  # Blank line to separate alignments
            cat "$ALN2"
        } > "$MERGE_ALN"
        echo "Generated merged file: $MERGE_ALN"

        # Copy control file into the output directory
        cp "$CONTROL_FILE" "$NEW_OUTPUT_DIR/"
        echo "Copied control file yn00.ctl to $NEW_OUTPUT_DIR"

        # Enter the new output directory
        pushd "$NEW_OUTPUT_DIR" > /dev/null
        echo "Changed working directory to $NEW_OUTPUT_DIR"

        # Run yn00 with the copied control file
        "$YN00_EXEC" "yn00.ctl" > yn00_output.txt
        echo "Finished running yn00; output saved to $NEW_OUTPUT_DIR/yn00_output.txt"

        # Return to the calling directory
        popd > /dev/null
        echo "----------------------------------------"

    done < "$COMPARISONS_FILE"
}

# Export the function for GNU parallel
export -f process_id

# Run GNU parallel with load/concurrency constraints
parallel --bar --jobs 20 --load 75% process_id {} ::: $(seq 1 1577)


##############
## 2. Summary ##
##############
#!/bin/bash
# Summaries Result files under each CDS_i directory (1-1577) using parallel workers.

# Process a single CDS directory
process_directory() {
    i="$1"

    # Paths for the current index
    OUTPUT_DIR="OUTPUT/CDS_${i}"
    SUMMARY_DIR="Summary/CDS_${i}/"

    # Ensure Summary directory exists
    mkdir -p "$SUMMARY_DIR"

    # Summary file paths
    SUMMARY1="${SUMMARY_DIR}/summary1.txt"
    SUMMARY2="${SUMMARY_DIR}/summary2.txt"
    SUMMARY3="${SUMMARY_DIR}/summary3.txt"
    SUMMARY4="${SUMMARY_DIR}/summary4.txt"

    # Initialize files with headers
    echo "(A) Nei & Gojobori 1986. dN/dS (dN, dS)" > "$SUMMARY1"
    echo "(B) Yang & Nielsen (2000) method" > "$SUMMARY2"
    echo "seq. seq.     S       N        t   kappa   omega     dN +- SE    dS +- SE" >> "$SUMMARY2"
    echo "(C) LWL85, LPB93 & LWLm methods" > "$SUMMARY3"
    echo "Number of codons" > "$SUMMARY4"

    # Walk through Result files
    find "$OUTPUT_DIR" -type f -name "Result" | while read -r result_file; do
        # Extract relevant sections
        first_output=$(grep -A 8 '(A) Nei-Gojobori' "$result_file" | tail -n 1)
        second_output=$(grep -A 9 '(A) Nei-Gojobori' "$result_file" | tail -n 1 | awk '{print $2,$3,$4}')
        third_output=$(grep -A 8 '(B) Yang & Nielsen' "$result_file" | tail -n 1 | cut -d " " -f 9-)
        forth_output=$(head -n 3 "$result_file" | tail -n 1)
        fifth_output=$(grep -A 16 '(C) LWL85, LPB93 & LWLm methods' "$result_file" | tail -n 8)

        # Append to summaries
        echo "${first_output} ${second_output}" >> "$SUMMARY1"
        echo "${first_output} ${third_output}" >> "$SUMMARY2"
        echo "${first_output} ${fifth_output}" >> "$SUMMARY3"
        echo "${first_output} ${forth_output}" >> "$SUMMARY4"
    done

    echo "Summaries for CDS_${i} written to the four files under Summary."
}

# Export to GNU parallel
export -f process_directory

# Run GNU parallel across all CDS directories
seq 1 1577 | parallel --jobs 5 process_directory {}
