#!/bin/bash

#! For each genome FASTA file of $strain

input_dir="../"


protein_db="../Hp26695.faa"


for fasta_file in "$input_dir"/*.fasta; do
    # Extract the file name (without path and extension)
    strain=$(basename "$fasta_file" .fasta)
    
    # Create the output directory named as strain'_prokka'
    output_dir="../Annotation/${strain}_prokka"
    mkdir -p "$output_dir"
    
    
    prokka --outdir "$output_dir" \
           --force \
           --rawproduct \
           --prefix "$strain" \
           --locustag "$strain" \
           --genus Helicobacter \
           --species pylori \
           --strain "$strain" \
           --kingdom Bacteria \
           --usegenus \
           --proteins "$protein_db" \
           --cpus 16 \
           --evalue 1e-24 \
           "$fasta_file"
    
    echo "Processing of $strain completed!"
done
