#!/usr/bin/env python3

"""
This scripts reads a nucleotide fasta file and concatenates the sequences into a single sequence.
Then it shuffles the sequence and writes it to a new fasta file.
"""

from Bio import SeqIO
from Bio.Seq import Seq
from random import shuffle
import sys
import os
import argparse
import re

# Get the parent directory of the script
script_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(script_dir)

# Add the parent directory to the system path
sys.path.append(parent_dir)
from fasta_get_ORFs import find_orfs

def concatenate_and_shuffle(input_file, min_size, num_orfs, strand='s', nocleanup=False):
    # Build the concatenated sequence once
    sequences = [str(record.seq) for record in SeqIO.parse(input_file, "fasta")]
    concatenated_sequence = ''.join(sequences)
    if not nocleanup:
        concatenated_sequence = re.sub('[^ATCG]', '', concatenated_sequence.upper())
    if strand == 'r':
        concatenated_sequence = str(Seq(concatenated_sequence).reverse_complement())

    orfs = []
    while len(orfs) < num_orfs:
        shuffled_sequence = list(concatenated_sequence)
        shuffle(shuffled_sequence)
        shuffled_sequence = ''.join(shuffled_sequence)
        if strand == 'r':
            shuffled_sequence = str(Seq(shuffled_sequence).reverse_complement())
        orfs += [orf for orf in find_orfs(shuffled_sequence, min_size, keep_first_orf=False, keep_last_orf=False) if len(orf[-1]) % 3 == 0]
    return orfs[:num_orfs]

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Concatenate and shuffle sequences from a fasta file.')
    parser.add_argument('input_file', type=str, help='Input fasta file')
    parser.add_argument('output_file', type=str, help='Output fasta name (NO EXTENSION)')
    parser.add_argument('--min_size', type=int, default=60, help='Minimum ORF size')
    parser.add_argument('--num_orfs', type=int, default=1, help='Number of ORFs to write to the output file')
    parser.add_argument('--strand', default='b', choices=['s', 'r', 'b'], help='Strand: forward (s), reverse (r), or both (b)')
    parser.add_argument('--nocleanup', action='store_true', help='Do not remove the unusual nucleotides')
    args = parser.parse_args()

    # Ensure min_size is set to 0 if provided as 0
    min_size = args.min_size if args.min_size != 0 else 0
 
    # If num_orfs is not provided, set it to the number of entries in the input fasta file
    if args.num_orfs is None:
        with open(args.input_file, 'r') as f:
            args.num_orfs = sum(1 for _ in SeqIO.parse(f, "fasta"))

    if args.strand == 'b':
        if args.num_orfs % 2 != 0:
            raise ValueError("The number of ORFs must be even when extracting from both strands")
        sense_num_orfs = args.num_orfs // 2
        antisense_num_orfs = args.num_orfs // 2
    else:
        sense_num_orfs = args.num_orfs
        antisense_num_orfs = args.num_orfs

    with open(args.output_file + ".fna", 'w') as f, open(args.output_file + ".faa", 'w') as f_translated:
        if args.strand in ['s', 'b']:
            sense_orfs = concatenate_and_shuffle(args.input_file, min_size, sense_num_orfs, strand='s', nocleanup=args.nocleanup)
            for i, orf in enumerate(sense_orfs):
                sequence = orf[-1]
                f.write(f'>orf_sense_{i+1}\n{sequence}\n')
                translated_sequence = str(Seq(sequence).translate())
                f_translated.write(f'>orf_sense_{i+1}\n{translated_sequence}\n')

        if args.strand in ['r', 'b']:
            antisense_orfs = concatenate_and_shuffle(args.input_file, min_size, antisense_num_orfs, strand='r', nocleanup=args.nocleanup)
            for i, orf in enumerate(antisense_orfs):
                sequence = orf[-1]
                f.write(f'>orf_antisense_{i+1}\n{sequence}\n')
                translated_sequence = str(Seq(sequence).translate())
                f_translated.write(f'>orf_antisense_{i+1}\n{translated_sequence}\n')