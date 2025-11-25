#!/bin/bash
set -e  # Exit immediately if any command fails
set -o pipefail  # Pipeline fails if any command within fails

# Determine directories relative to this script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

# Define output directory
OUTPUTDIR="${BASE_DIR}/OUTPUT/WGS"
cd "${OUTPUTDIR}" || { echo "Unable to switch to directory ${OUTPUTDIR}"; exit 1; }

# Record start time
date > WGS_HP.time
echo "Analysis started at $(date)" >> WGS_HP.time

# Run the initial fineSTRUCTURE command
echo "Running initial fineSTRUCTURE command..."
fs WGS_HP.cp -hpc 1 \
    -idfile WGS_HP.ids \
    -phasefiles WGS_HP.phase \
    -recombfiles WGS_HP.recombfile \
    -s3iters 200000 \
    -s4iters 50000 \
    -s1minsnps 1000 \
    -s1indfrac 0.1 \
    -go

# Define helper to execute command files
run_commandfile() {
    local cmdfile=$1
    local stage=$2

    if [ -f "${cmdfile}" ]; then
        echo "Running Stage ${stage} command file: ${cmdfile}"
        parallel -j 16 < "${cmdfile}" || { echo "Parallel execution failed for Stage ${stage}."; exit 1; }
        echo "Stage ${stage} completed. Continuing analysis..."
        fs WGS_HP.cp -go || { echo "fs command failed after Stage ${stage}."; exit 1; }
    else
        echo "Warning: ${cmdfile} does not exist. Skipping Stage ${stage}."
    fi
}

# Run Stage 1 through Stage 4
for stage in {1..4}; do
    cmdfile="WGS_HP/commandfiles/commandfile${stage}.txt"
    run_commandfile "${cmdfile}" "${stage}"
done

# Record end time
date >> WGS_HP.time
echo "Analysis ended at $(date)" >> WGS_HP.time

echo "fineSTRUCTURE analysis completed."
