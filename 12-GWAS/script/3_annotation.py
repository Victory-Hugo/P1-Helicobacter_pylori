#!/usr/bin/env python3
import os
import argparse
import pandas as pd
import numpy as np

def extract_product(attr_str):
    parts = attr_str.split("product=")
    return parts[1].split(";")[0] if len(parts) > 1 else None

def gwas_analysis(gwas_file, gff_file, output_dir, logp_threshold=5, dist_threshold=1000,
                  gwas_sep="\t", gff_sep="\t", gwas_header=0, gff_header=None):
    """
    Load GWAS and GFF files, process the data, and write annotated output files.

    :param gwas_file: Path to the GWAS result file
    :param gff_file:  Path to the GFF annotation file
    :param output_dir: Output directory
    :param logp_threshold: -log10(p-value) threshold used to keep significant SNPs
    :param dist_threshold: Distance threshold (bp) to consider SNPs as being close
    :param gwas_sep: Delimiter for the GWAS file
    :param gff_sep: Delimiter for the GFF file
    :param gwas_header: Header argument passed to pandas for the GWAS file
    :param gff_header: Header argument passed to pandas for the GFF file (use None if no header)
    """
    # Create the output directory
    os.makedirs(output_dir, exist_ok=True)

    # 1. Load GWAS data
    gwas_data = pd.read_csv(gwas_file, sep=gwas_sep, header=gwas_header)
    print("Columns in GWAS data:", gwas_data.columns.tolist())

    # 2. Handle -log10(p-value) information
    if 'negLog10' in gwas_data.columns:
        gwas_data['logp'] = gwas_data['negLog10']
    elif 'pvalue' in gwas_data.columns:
        gwas_data['logp'] = -np.log10(gwas_data['pvalue'])
    else:
        raise ValueError("The input file must contain either a 'negLog10' or 'pvalue' column.")

    # 3. Confirm that column 'ps' exists for SNP positions
    if 'ps' not in gwas_data.columns:
        raise ValueError("Column 'ps' is required to determine SNP positions.")

    # 4. Filter significant SNPs with finite -log10(p-value) scores above the threshold
    gwas_data.replace([np.inf, -np.inf], np.nan, inplace=True)
    sig = gwas_data.loc[(gwas_data['logp'] > logp_threshold) & (gwas_data['logp'].notna())].copy()

    # 5. Sort by SNP position, compute distances, and keep entries closer than the threshold
    sig.sort_values('ps', inplace=True)
    sig.reset_index(drop=True, inplace=True)
    sig['delta'] = sig['ps'].diff().abs().fillna(0)
    sig = sig.loc[sig['delta'] < dist_threshold].drop(columns='delta')

    # 6. Load CDS annotations from the GFF file
    gff_CDS = pd.read_csv(gff_file, sep=gff_sep, header=gff_header,
                            names=['seqid', 'source', 'type', 'start', 'end', 'score',
                                   'strand', 'phase', 'attributes'])

    # 7. Match each significant SNP to CDS annotations and record the negLog10 value
    sig_CDS_list = []
    for idx, snp_row in sig.iterrows():
        pos = snp_row['ps']
        current_logp = snp_row['logp']  # Track the current SNP -log10(p-value)
        matched = gff_CDS.loc[(gff_CDS['start'] < pos) & (gff_CDS['end'] > pos)]
        if not matched.empty:
            matched = matched.copy()
            matched['negLog10'] = current_logp
            sig_CDS_list.append(matched)
    if sig_CDS_list:
        # Combine all matches
        sigCDS_all = pd.concat(sig_CDS_list, ignore_index=True)
    else:
        sigCDS_all = pd.DataFrame(columns=list(gff_CDS.columns) + ['negLog10'])

    # 8. Extract product details from the attributes field
    sigCDS_all['product'] = sigCDS_all['attributes'].apply(extract_product)

    # 9. Keep the maximum negLog10 per CDS (defined by start, end, product)
    sigCDS_grouped = sigCDS_all.groupby(['start', 'end', 'product'], as_index=False)['negLog10'].max()

    # Sort by CDS start position and export start, end, product, negLog10
    sigCDS_sorted = sigCDS_grouped.sort_values('start')
    cds_output_path = os.path.join(output_dir, "significantCDS.csv")
    sigCDS_sorted.to_csv(cds_output_path, index=False, header=False, sep=",")

    # 10. Write detailed SNP information
    snp_output_path = os.path.join(output_dir, "significantSNPs.txt")
    sig.to_csv(snp_output_path, sep='\t', index=False)

    print("Analysis complete. Results written to:", output_dir)

def main():
    parser = argparse.ArgumentParser(description="GWAS annotation helper that can be configured from shell variables")
    parser.add_argument("--gwas_file", type=str, required=True, help="Path to the GWAS result file")
    parser.add_argument("--gff_file", type=str, required=True, help="Path to the CDS annotation GFF file")
    parser.add_argument("--output_dir", type=str, required=True, help="Output directory")
    parser.add_argument("--logp_threshold", type=float, default=5, help="-log10(p-value) filtering threshold (default: 5)")
    parser.add_argument("--dist_threshold", type=int, default=1000, help="Proximity threshold between SNPs in bp (default: 1000)")
    parser.add_argument("--gwas_sep", type=str, default="\t", help="GWAS file delimiter (default: tab)")
    parser.add_argument("--gff_sep", type=str, default="\t", help="GFF file delimiter (default: tab)")
    parser.add_argument("--gwas_header", type=int, default=0, help="Header argument for the GWAS file (default: 0)")
    parser.add_argument("--gff_header", type=str, default="None", help="Header argument for the GFF file; use None if there is no header")
    args = parser.parse_args()

    # Convert gff_header argument when strings are used
    gff_header = None if args.gff_header == "None" else int(args.gff_header)

    gwas_analysis(
        gwas_file=args.gwas_file,
        gff_file=args.gff_file,
        output_dir=args.output_dir,
        logp_threshold=args.logp_threshold,
        dist_threshold=args.dist_threshold,
        gwas_sep=args.gwas_sep,
        gff_sep=args.gff_sep,
        gwas_header=args.gwas_header,
        gff_header=gff_header
    )

if __name__ == "__main__":
    main()
