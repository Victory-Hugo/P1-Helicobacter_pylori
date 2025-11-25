#!/bin/bash



input_dir="../Annotation"


output_dir="../Alignment/"

# 创建输出目录
mkdir -p "$output_dir"

# Collect all .gff files
input_files=$(find "$input_dir" -type f -name "*.gff")

# Exit if no .gff files are found
if [ -z "$input_files" ]; then
    echo "No .gff files found in the specified directory!"
    exit 1
fi


panaroo -i $input_files \
        -o "$output_dir/panaroo_output/" \
        --clean-mode strict \
        --threshold 0.90 \
        -a core \
        --core_threshold 0.95 \
        --len_dif_percent 0.75 \
        -t 20

echo "Panaroo analysis completed!"
