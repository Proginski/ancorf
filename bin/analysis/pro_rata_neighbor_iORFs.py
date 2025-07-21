#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This program reads a fasta where headers are genome phylip names.
It also reads a two columns file where the first column is the genome name and the second column the phylip name.
Then if neighbor_A is found in 12% of the headers of the FASTA, it picks 120 random sequences from **args.iorf_dir**/**neighbor_A**_iORF_1000.fna .
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
    parser.add_argument('--name_mapping', required=True, help='The two columns file with genome names and phylip names.')
    parser.add_argument('--iorf_dir', required=True, help='The directory with the iORFs FASTA files.')
    args = parser.parse_args()

    # Read the name_mapping file into a pandas DataFrame
    names_df = pd.read_csv(args.name_mapping, sep='\t', header=None, names=['genome_name', 'phylip_name'])

    # Read the input_fasta file and count the number of sequences for each phylip_name
    sequence_counts = {}
    for record in SeqIO.parse(args.input_fasta, "fasta"):
        phylip_name = record.id.split("__")[0]
        if phylip_name in sequence_counts:
            sequence_counts[phylip_name] += 1
        else:
            sequence_counts[phylip_name] = 1

    # Distribute the 1000 output sequences proportionally based on the counts
    # Normalize the counts to ensure they sum to 1000
    total_sequences = sum(sequence_counts.values())
    sequence_counts = {k: round(v / total_sequences * 1000) for k, v in sequence_counts.items()}
    total_sequences = sum(sequence_counts.values())
    if total_sequences < 1000:
        num_sequences_to_add = 1000 - total_sequences
        for i in range(num_sequences_to_add):
            # Pick a random key in sequence_counts to add an extra sequence
            random_key = random.choice(list(sequence_counts.keys()))
            sequence_counts[random_key] += 1
    elif total_sequences > 1000:
        num_sequences_to_remove = total_sequences - 1000
        for i in range(num_sequences_to_remove):
            # Pick a random key in sequence_counts to remove a sequence
            volonteer = False
            while not volonteer:
                random_key = random.choice(list(sequence_counts.keys()))
                if sequence_counts[random_key] > 0:
                    volonteer = True
                    sequence_counts[random_key] -= 1
    
    # For each phylip_name, pick a number of random sequences from the corresponding iORF file
    total_count = sum(sequence_counts.values())
    output_records = []
    for index, row in names_df.iterrows():
        genome_name = row['genome_name']
        phylip_name = row['phylip_name']
        if phylip_name in sequence_counts:
            num_to_select = sequence_counts[phylip_name]
            iorf_file = f"{args.iorf_dir}/{genome_name}_iORF_1000.fna"
            iorf_records = list(SeqIO.parse(iorf_file, "fasta"))
            if len(iorf_records) >= num_to_select:
                random_records = random.sample(iorf_records, num_to_select)
            else:
                random_records = iorf_records
            for record in random_records:
                record.id = f"{record.id}__{phylip_name}"
            output_records.extend(random_records)

    # Write the output_records to stdout
    SeqIO.write(output_records, sys.stdout, "fasta")

if __name__ == '__main__':
    main()