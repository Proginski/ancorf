#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import glob
import argparse
from Bio import SeqIO

def sanitize_filename(filename):
    return filename.replace(":", "__COLON__").replace("/", "__SLASH__").replace("!", "__EXCLAMATION__").replace("|", "__PIPE__")

def extract_sequences(input_dir, output_dir, mode, evalue_threshold=None, aligner="lalign36"):
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Find all TSV files in the input directory
    tsv_files = glob.glob(f"{input_dir}/*_{aligner}.tsv")

    # Initialize dictionaries to store sequences
    nucl_sequences = {}
    aa_sequences = {}

    # Process each TSV file
    for tsv_file in tsv_files:
        found=False
        with open(tsv_file, 'r') as f:
            for line in f:
                # Skip comment lines
                if line.startswith('#'):
                    continue

                # Extract query name, target and e-value from the line
                query, target, _, _, _, _, _, _, _, _, evalue, _ = line.split('\t')
                ORFid = f"{query}__ANC_FRAG__{target}"
                evalue = float(evalue)

                # Skip the hits that do not meet the e-value threshold
                if evalue > evalue_threshold:
                    continue

                # Find the corresponding FASTA files
                sanitized_query = sanitize_filename(query)
                nucl_fasta_file = f"{input_dir}/{sanitized_query}_anc_ORFs.fna"
                aa_fasta_file = f"{input_dir}/{sanitized_query}_anc_ORFs.faa"

                # Extract the sequence with the matching ID from the FASTA files
                for record in SeqIO.parse(nucl_fasta_file, "fasta"):
                    if record.id == target:
                        nucl_sequences[ORFid] = record
                        nucl_sequences[ORFid].id = ORFid
                        found=True
                        break
                for record in SeqIO.parse(aa_fasta_file, "fasta"):
                    if record.id == target:
                        aa_sequences[ORFid] = record
                        aa_sequences[ORFid].id = ORFid
                        break

                # If mode is 'best', only process the best (first) non-comment line
                if mode == 'best':
                    break

            if not found:
                print(f"No hits found for {tsv_file}")

    # Write the output FASTA files
    SeqIO.write(nucl_sequences.values(), f"{output_dir}/anc_ORFs.fna", "fasta")
    SeqIO.write(aa_sequences.values(), f"{output_dir}/anc_ORFs.faa", "fasta")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Extract sequences from FASTA files.')
    parser.add_argument('input_dir', help='The directory containing the input TSV and FASTA files.')
    parser.add_argument('output_dir', help='The directory where the output FASTA files will be written.')
    parser.add_argument('--mode', choices=['best', 'threshold'], default='best', help='The mode of operation.')
    parser.add_argument('--evalue', type=float, default=1e-2, help='The e-value threshold for selecting hits (only used in threshold mode).')
    parser.add_argument('--aligner', choices=['lalign36', 'ssearch36'], default='lalign36', help='The aligner used to generate the TSV files.')
    args = parser.parse_args()

    # Call the function with the input and output directories, mode and e-value threshold provided as command-line arguments
    extract_sequences(args.input_dir, args.output_dir, args.mode, args.evalue, args.aligner)