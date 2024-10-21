#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This program takes a nucleotide FASTA and extracts all its ORFs from STOP codon to STOP codon.
It includes the last STOP codon.
"""

import argparse
import sys
from Bio import SeqIO
from Bio.Seq import Seq

def find_orfs(seq, min_orf_length, keep_first_orf=True, keep_last_orf=True):

    if min_orf_length < 0:
        raise ValueError("Minimum protein length must be a positive integer")
    
    seq = seq.upper()  # Convert the sequence to uppercase
    stop_codons = ["TAA", "TAG", "TGA"]

    orf_list = []

    for frame in range(3):
        start = frame
        for i in range(frame, len(seq), 3):

            ORF_end_i = "NA" # If ORF_end_i is "NA", it means that the ORF is not complete

            current_last_codon = seq[i:i+3]

            if current_last_codon in stop_codons: # if the current last codon is a STOP
                ORF_end_i = i+3
            elif len(current_last_codon) != 3: # if the current last codon is not complete (end of the sequence)
                ORF_end_i = i

            if ORF_end_i != "NA": # if the ORF is complete
                if keep_first_orf or start != frame:
                    orf = seq[start:ORF_end_i]
                    if len(orf) - 3 >= min_orf_length:  # Subtract 3 for the stop codon to check the ORF length
                        if current_last_codon in stop_codons or keep_last_orf:
                            orf_list.append((frame, start, ORF_end_i, orf))
        
            if current_last_codon in stop_codons: # if the current last codon is a STOP
                start = i + 3
            elif len(current_last_codon) != 3: # if the current last codon is not complete (end of the sequence)
                break

    return orf_list

def main():
    parser = argparse.ArgumentParser(description='Extract ORFs from FASTA')
    parser.add_argument('input', help='Input FASTA file')
    parser.add_argument('-o', '--output', help='Output name (optional, without extension)')
    parser.add_argument('-m', '--min_size', type=int, default=60, help='Minimum ORF size (STOP excluded, default: 60)')
    parser.add_argument('-t', '--type', default='both', choices=['nucl', 'aa', 'both'], help='Output type: nucleotides (nucl), their translation (aa), or both')
    parser.add_argument('-s', '--strand', default='b', choices=['s', 'r', 'b'], help='Strand: forward (s), reverse (r), or both (b)')

    args = parser.parse_args()

    fasta_sequences = SeqIO.parse(open(args.input),'fasta')

    output_nucl = open(f"{args.output}.fna", "w") if args.output else sys.stdout
    output_aa = open(f"{args.output}.faa", "w") if args.output else sys.stdout

    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        if args.strand in ['s', 'b']:
            orfs = find_orfs(sequence, args.min_size)
            for frame, start, end, orf in orfs:
                if args.type in ["nucl", "both"]:
                    print(f'>{name}_frame{frame+1}_{start+1}-{end}\n{orf}', file=output_nucl)
                if args.type in ["aa", "both"]:
                    print(f'>{name}_frame{frame+1}_{start+1}-{end}\n{Seq(orf).translate()}', file=output_aa)

        if args.strand in ['r', 'b']:
            sequence = str(fasta.seq.reverse_complement())
            orfs = find_orfs(sequence, args.min_size)
            for frame, start, end, orf in orfs:
                if args.type in ["nucl", "both"]:
                    print(f'>{name}_frame{frame+1}_{start+1}-{end}_reverse\n{orf}', file=output_nucl)
                if args.type in ["aa", "both"]:
                    print(f'>{name}_frame{frame+1}_{start+1}-{end}_reverse\n{Seq(orf).translate()}', file=output_aa)

    if args.output:
        output_nucl.close()
        output_aa.close()

if __name__ == "__main__":
    main()