# 1-nucmer: Helicobacter pylori Variant Calling Pipeline

This directory contains a set of configuration files and scripts to run a whole-genome alignment and SNP-calling pipeline for *Helicobacter pylori* using `nucmer` (from MUMmer) and downstream processing tools.

The pipeline takes genome assemblies (or contigs) as input, aligns them to a reference genome, extracts SNPs, and produces merged VCF files and core-genome SNP sets for population or phylogenetic analysis.

## Directory Structure

- `conf/`
  - `NC_000915.fasta`: Reference genome sequence for *H. pylori* (e.g. strain 26695).
  - `NC_000915.csv`: Annotation or coordinate information linked to the reference genome.
  - `N_process/`: Configuration and helper data for handling Ns (masked or low-quality bases) during processing.
- `script/`: Shell and Python scripts implementing each step of the pipeline.

## Main Workflow Scripts

All scripts below are located in the `script/` directory and are intended to be run in sequence (or via a wrapper script, if provided).

- `1-nucmer.sh`  
  Runs `nucmer` to align each input genome/assembly to the reference (`NC_000915.fasta`), and generates alignment files (e.g. `.delta` / `.coords`) for downstream SNP calling.

- `2_tsv_df.py`  
  Parses the raw `nucmer`/`show-snps` output (typically TSV-like text) and converts it into a cleaned, tabular data frame (TSV/CSV), standardizing columns such as position, reference allele, and alternative allele.

- `3_coord_csv.py`  
  Converts coordinate-based outputs (e.g. from `show-coords` or similar) into CSV format, linking alignment coordinates to reference positions and facilitating merging with SNP tables.

- `3_fillN.py`  
  Handles positions where the input genome contains `N` or low-quality bases. This script can mark, filter, or replace such positions to avoid false SNP calls.

- `3_trim_csv.py`  
  Trims or filters the CSV tables (e.g. removing low-coverage sites, overlapping regions, or non-biallelic variants) to produce high-confidence SNP tables.

- `4_vcf.py`  
  Converts the processed SNP tables into VCF format for each sample, preserving reference coordinates and variant annotations.

- `5-1_merge.sh`  
  Merges per-sample VCF or SNP tables into a combined matrix (multi-sample VCF or SNP table) across all genomes.

- `5-2_merge.sh`  
  Performs additional merging or post-processing (e.g. merging multiple batches, harmonizing sample IDs, or reordering sites) after the initial merge.

- `6_0_N_processC++.sh`  
  Wrapper script to run C++-based tools for more efficient `N`-processing or complex filtering steps on large SNP matrices.

- `7-extractSNP.sh`  
  Extracts SNP-only positions from the merged data (e.g. removing indels and non-variable sites) to obtain a clean SNP alignment or SNP list for downstream analyses.

- `8-core_genome.sh`  
  Identifies core-genome SNPs shared across all (or a defined majority of) samples and outputs a core-genome SNP alignment suitable for phylogenetic tree construction or population structure analysis.

- `run_filter.sh`  
  A high-level driver script that chains several filtering steps together (e.g. trimming, N-handling, SNP extraction) according to predefined criteria.

## Typical Usage

1. Prepare input genomes and place them in the expected input directory (as defined inside the scripts).
2. Ensure `nucmer` (MUMmer), Python, and any required C++ tools or libraries are installed and accessible in `PATH`.
3. Edit configuration paths in the scripts or `conf/` files if necessary (e.g. reference file locations, input/output directories).
4. Run the pipeline step by step:
   - `bash script/1-nucmer.sh`
   - `python script/2_tsv_df.py`
   - `python script/3_coord_csv.py`
   - `python script/3_fillN.py`
   - `python script/3_trim_csv.py`
   - `python script/4_vcf.py`
   - `bash script/5-1_merge.sh`
   - `bash script/5-2_merge.sh`
   - `bash script/6_0_N_processC++.sh`
   - `bash script/7-extractSNP.sh`
   - `bash script/8-core_genome.sh`
5. Use `run_filter.sh` to reapply filtering or to customize filtering criteria as needed.

## Output

The final outputs typically include:

- Per-sample VCF files aligned to `NC_000915.fasta`.
- A merged multi-sample SNP matrix or VCF.
- SNP-only and core-genome SNP alignments suitable for phylogenetic and population genetic analyses of *Helicobacter pylori*.

