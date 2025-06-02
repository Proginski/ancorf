#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script takes as input a nucleotide FASTA file and outputs two fasta files :
    one with the shuffled input sequences with no stop but the original ending stop
    a second, which is the AA translation of the first one.
    
To conceal run time and codon pseudo-randomness, every generated stop codon is 
shuffled together with the n flanking codons (default : 1 at each side).
"""

import random
import sys
from copy import deepcopy
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse

def check_fasta_file(fasta_file):
    with open(fasta_file, 'r') as file:
        # Read the file
        content = file.read()
        # Check if there are any entries
        if '>' not in content:
            print(f"Error: {fasta_file} has no entries", file=sys.stderr)
            sys.exit(1)

def generate_unique_ids(records):
    """Assign unique IDs to all records without modifying the original objects."""
    occ_dict = {}
    new_records = []
    for record in records:
        # Check if the record ID already exists in the dictionary
        if record.id in occ_dict:
            occ_dict[record.id] += 1
        else:
            occ_dict[record.id] = 1
        # Create a deep copy of the record and assign a new ID
        new_record = deepcopy(record)
        new_record.id = f"{record.id}_{occ_dict[record.id]}"
        new_record.description = ""
        new_records.append(new_record)
    return new_records

def shuffle_out(seq, side_codons=1):
    stop_codons = ['TAA', 'TAG', 'TGA']
    #print(f"Initial sequence: {seq} and its translation: {Seq(seq).translate()}")
    seq = ''.join(random.sample(seq, len(seq))) # Initial shuffling
    stop_index = -1
    for i in range(0, len(seq) - 2, 3):
        if seq[i:i+3] in stop_codons:
            stop_index = i
            break
    while stop_index != -1:
        # If the sequence has less than one stop codon plus two time the number of neighbor codon (because n on the left, n on the ritgh), we shuffle the whole sequence
        if len(seq) <= 3*(side_codons*2 +1) :
            seq = ''.join(random.sample(seq, len(seq)))
            #print("Shuffling the whole sequence")
        else :
            # First try to get n codons on the left and n on the right
            if stop_index - 3*side_codons >= 0 and stop_index + 3*(1+side_codons) <= len(seq):
                segment_start = stop_index - 3*side_codons
            # Else try to get 2n codons on the right
            elif stop_index + 3*(1+side_codons) <= len(seq):
                segment_start = stop_index
            # Else try to get 2n codons on the left
            elif stop_index - 3*(2*side_codons) >= 0:
                segment_start = stop_index - 3*(2*side_codons)
            segment = seq[segment_start:segment_start + 3*(1+2*side_codons)]
            
            #print(f"Shuffling segment {segment} (translation: {Seq(segment).translate()})")
            segment = ''.join(random.sample(segment, len(segment)))
            seq = seq[:segment_start] + segment + seq[segment_start+len(segment):]
            
        stop_index = -1
        #print("Intermediate sequence: ", seq , " and its translation: ", Seq(seq).translate())
        for i in range(0, len(seq) - 2, 3):
            if seq[i:i+3] in stop_codons:
                stop_index = i
                break
    #print(f"Final sequence: {seq} and its translation: {Seq(seq).translate()}")
    return seq

def shuffle_sequence(input_file, stop, nb_seq, side_codons, output):
    # Read the sequences from the input file
    records = list(SeqIO.parse(input_file, "fasta"))
    # If nb_seq is specified, generate the desired number of sequences
    if nb_seq is not None:
        if nb_seq == len(records):
            pass
        elif nb_seq < len(records):
            records = random.sample(records, nb_seq)
        elif nb_seq > len(records):
            while len(records) < nb_seq:
                records.extend(random.sample(records, min(len(records), nb_seq - len(records))))

    # Assign unique IDs
    records = generate_unique_ids(records)

    # Convert all sequences to uppercase and shuffle them
    for record in records:
        record.seq = Seq(str(record.seq).upper())
        if stop == 'in':
            # Just shuffle the sequence (there may be stop codons in it)
            record.seq = Seq(''.join(random.sample(str(record.seq), len(record.seq))))
        else:
            # Store the terminal stop codon and shuffle the rest of the sequence avoiding inframe stops
            # Check if the sequence ends with a stop codon
            if str(record.seq[-3:]) in ['TAA', 'TAG', 'TGA']:
                stop_codon = str(record.seq[-3:])
                rest_of_seq = str(record.seq[:-3])
            else:
                stop_codon = ''
                rest_of_seq = str(record.seq)
            shuffled_seq = shuffle_out(rest_of_seq, side_codons)
            record.seq = Seq(shuffled_seq + stop_codon)
            record.description = ""
            record.id = f"{record.id}_random_{stop}"

    # Write the shuffled sequences to the output file
    SeqIO.write(records, f"{output}.fna", "fasta")

    # Convert the sequences to amino acids and write to the output file
    aa_records = []
    for record in records:
        aa_records.append(SeqRecord(record.seq.translate(), id=record.id, description=record.description))
    # aa_records = [SeqRecord(record.seq.translate(), id=record.id, description=record.description) for record in records]
    SeqIO.write(aa_records, f"{output}.faa", "fasta")

def main():
    parser = argparse.ArgumentParser(description='Shuffle sequences in a FASTA file.')
    parser.add_argument('input_file', type=str, help='Input FASTA file.')
    parser.add_argument('--stop', type=str, choices=['in', 'out'], default='out',
                        help='Specify if this is a simple shuffling (in) or if we avoid to generate inframe stops (out).')
    parser.add_argument('--side_codons', type=int, default=1, help='Specify the number of codons to shuffle on each side of a stop codon (only used with "--stop out").')
    parser.add_argument('--output', type=str, default='output', help='Specify the base name for the output FASTA files.')
    parser.add_argument('--nb_seq', type=int, default=None, help='Specify the number of new sequences. Default is as many as the input FASTA.')

    args = parser.parse_args()

    # Check if the input file is a valid FASTA file
    check_fasta_file(args.input_file)

    # If nb_seq is not specified, read the sequences from the input file to determine the default value for nb_seq
    if args.nb_seq is None:
        records = list(SeqIO.parse(args.input_file, "fasta"))
        args.nb_seq = len(records)

    print(f"Shuffling {args.nb_seq} sequences from {args.input_file} to {args.output}.fna and {args.output}.faa")
    shuffle_sequence(args.input_file, args.stop, args.nb_seq, args.side_codons, args.output)

if __name__ == "__main__":
    main()