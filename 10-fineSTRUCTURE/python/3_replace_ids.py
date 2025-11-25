#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import csv
import argparse

def load_id_anchor_map(csv_file):
    """
    Load the ID→Anchor mapping from a CSV file.
    Expectation: the first row is a header, ID is column 1, Anchor is column 2, comma-separated.
    """
    id_anchor_dict = {}
    with open(csv_file, 'r', encoding='utf-8') as f:
        reader = csv.reader(f)
        header = next(reader, None)  # Skip the header
        for row in reader:
            if not row:
                continue
            the_id = row[0].strip()
            anchor = row[1].strip()
            id_anchor_dict[the_id] = anchor
    return id_anchor_dict

def process_ids_file(ids_file, id_anchor_dict):
    """
    Read the .ids file line by line, look up the Anchor for column 1,
    and replace column 2 with the Anchor; keep the original if no match.
    """
    temp_file = ids_file + ".tmp"
    with open(ids_file, 'r', encoding='utf-8') as fin, \
         open(temp_file, 'w', encoding='utf-8') as fout:
        
        for line in fin:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            # Expected format: column 1 = ID, column 2 = default (usually same as ID), column 3 = additional info
            if len(parts) < 3:
                # Skip malformed lines
                continue
            
            the_id = parts[0]
            original_val = parts[1]
            col3 = parts[2]
            
            if the_id in id_anchor_dict:
                fout.write(f"{the_id} {id_anchor_dict[the_id]} {col3}\n")
            else:
                fout.write(f"{the_id} {original_val} {col3}\n")
    
    # Replace the original file with the processed temp file
    os.replace(temp_file, ids_file)

def main():
    parser = argparse.ArgumentParser(
        description="Process .ids files with a CSV mapping so the second column becomes the corresponding Anchor value"
    )
    parser.add_argument("--inf_csv", required=True,
                        help="Path to CSV containing ID and Anchor mappings")
    parser.add_argument("ids_files", nargs="+",
                        help="One or more .ids file paths to process")
    args = parser.parse_args()

    id_anchor_dict = load_id_anchor_map(args.inf_csv)
    
    for ids_file in args.ids_files:
        if os.path.isfile(ids_file):
            process_ids_file(ids_file, id_anchor_dict)
        else:
            print(f"Warning: {ids_file} is not a valid file")

if __name__ == "__main__":
    main()
