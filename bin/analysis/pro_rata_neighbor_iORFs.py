#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This program reads a fasta where headers are genome phylip names.
It also reads a two columns file where the first column is the genome name and the second column the phylip name.
Then if neighbor_A is found in 12% of the headers of the FASTA, it picks 120 random sequences from ../IORF/neighbor_A_iORF_1000.fna .
At the end it print a new FASTA with 1000 entries.
"""

import argparse
import random
from Bio import SeqIO
import pandas as pd
import sys

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_fasta', help='The FASTA file with genome phylip names in headers.')
    parser.add_argument('names_to_phylip_names', help='The two columns file with genome names and phylip names.')
    args = parser.parse_args()

    # Read the names_to_phylip_names file into a pandas DataFrame
    names_df = pd.read_csv(args.names_to_phylip_names, sep='\t', header=None, names=['genome_name', 'phylip_name'])

    # Read the input_fasta file and count the number of sequences for each phylip_name
    sequence_counts = {}
    for record in SeqIO.parse(args.input_fasta, "fasta"):
        phylip_name = record.id
        if phylip_name in sequence_counts:
            sequence_counts[phylip_name] += 1
        else:
            sequence_counts[phylip_name] = 1

    # For each phylip_name, pick a number of random sequences from the corresponding iORF file
    total_count = sum(sequence_counts.values())
    output_records = []
    for index, row in names_df.iterrows():
        genome_name = row['genome_name']
        phylip_name = row['phylip_name']
        if phylip_name in sequence_counts:
            count = sequence_counts[phylip_name]
            proportion = count / total_count
            num_to_select = round(proportion * 1000)
            iorf_file = f"../IORF/{genome_name}_iORF_1000.fna"
            iorf_records = list(SeqIO.parse(iorf_file, "fasta"))
            if len(iorf_records) >= num_to_select:
                random_records = random.sample(iorf_records, num_to_select)
            else:
                random_records = iorf_records
            output_records.extend(random_records)

    # Write the output_records to stdout
    SeqIO.write(output_records, sys.stdout, "fasta")

if __name__ == '__main__':
    main()