#!/bin/bash
# Generate 16 TreeMix run scripts for migration edges m=0..15, with 10 replicates per m
#*===============================
# Determine project directories relative to this script
SCRIPT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_ROOT}/.." && pwd)"

# todo Set the input file path; generated via ./python/plink2treemix.py
INPUT_FILE="${PROJECT_ROOT}/data/TreeMix_China3219.gz"
# todo Set the output directory
OUTPUT_DIR="${PROJECT_ROOT}/output/China_only"
# todo Set the directory to store generated scripts
SCRIPT_DIR="${PROJECT_ROOT}/script/China_only"
# todo Set the outgroup
OUTGROUP="hpAfrica2"
#*===============================
echo "Generating TreeMix analysis scripts..."

# Ensure output directories exist
mkdir -p "$OUTPUT_DIR"
mkdir -p "$SCRIPT_DIR"



# Generate scripts for each migration edge count
for m in {0..15}; do
  script_name="${SCRIPT_DIR}/treemix_m${m}.sh"
  
  cat > "$script_name" << EOF
#!/bin/bash
#
# TreeMix analysis script - migration edge count m=${m}, replicates rep=1..10
# Usage: bash treemix_m${m}.sh

set -euo pipefail

INPUT_FILE="$INPUT_FILE"
OUTPUT_DIR="$OUTPUT_DIR"
OUTGROUP="$OUTGROUP"

echo "Starting analysis for migration edges m=${m}..."

mkdir -p "\$OUTPUT_DIR"

# Perform 10 replicates per migration edge count
for rep in {1..10}; do
  prefix="\${OUTPUT_DIR}/Treemix${m}.\${rep}"
  echo "  [m=${m} rep=\${rep}] Running TreeMix..."
  treemix -i "\$INPUT_FILE" \\
    -root "\$OUTGROUP" \\
    -o "\$prefix" \\
    -m ${m} \\
    -se -bootstrap \\
    -global \\
    -k 500 \\
    -threads 8

  echo "  [m=${m} rep=\${rep}] Analysis completed"
done

echo "All replicates completed: m=${m}"
EOF

  chmod +x "$script_name"
  echo "Generated script: $script_name"
done

echo "All scripts generated in $SCRIPT_DIR"
