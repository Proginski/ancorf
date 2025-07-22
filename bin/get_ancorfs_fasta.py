#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This script builds FASTA files with ancestral ORFs.
get_ancorfs_fasta.py $ssearch36_tsvs --mode best --mode 1e-3 --mode 1e-2
Here are the fields of the TSV files:
# Fields: query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
If mode is 'best', the script will keep the best hit (best evalue).
If mode is a float, the script will keep the hits with an evalue <= mode.
"""

import argparse
import csv
import os

def sanitize_filename(filename):
    return filename.replace(":", "__COLON__").replace("/", "__SLASH__").replace("!", "__EXCLAMATION__").replace("|", "__PIPE__").replace("(", "__OPEN_PAREN__").replace(")", "__CLOSE_PAREN__").replace("+", "__PLUS__")

def read_fasta(fasta_file):
    sequences = {}
    try:
        with open(fasta_file, 'r') as file:
            sequence_id = None
            sequence = []
            for line in file:
                line = line.strip()
                if line.startswith('>'):
                    if sequence_id is not None:
                        sequences[sequence_id] = ''.join(sequence)
                    sequence_id = line[1:]
                    sequence = []
                else:
                    sequence.append(line)
            if sequence_id is not None:
                sequences[sequence_id] = ''.join(sequence)
    except FileNotFoundError:
        print(f"FASTA file {fasta_file} not found.")
    return sequences

def process_tsv_file(tsv, mode, fna_file, faa_file, written_pairs):
    with open(tsv, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        
        fieldnames = [
            'query_id', 'subject_id', 'identity', 'alignment_length', 
            'mismatches', 'gap_opens', 'q_start', 'q_end', 's_start', 
            's_end', 'evalue', 'bit_score'
        ]
        
        best_hit = None

        for row in reader:
            # Skip comment rows
            if row[0].startswith('#'):
                continue
            
            # Check if the row has the expected number of fields
            if len(row) != len(fieldnames):
                print(f"Skipping malformed row: {row}")
                continue

            row_dict = dict(zip(fieldnames, row))
            evalue = row_dict['evalue']
            if evalue is None:
                continue
            evalue = float(evalue)
            pair = (row_dict['query_id'], row_dict['subject_id'])
            if mode == 'best':
                if best_hit is None or evalue < float(best_hit['evalue']):
                    best_hit = row_dict
            else:
                if evalue <= float(mode) and pair not in written_pairs:
                    query_id = row_dict['query_id']
                    fasta_sequences = read_fasta(f"{sanitize_filename(query_id)}_anc_ORFs.fna")
                    translated_sequences = read_fasta(f"{sanitize_filename(query_id)}_anc_ORFs.faa")
                    if row_dict['subject_id'] in fasta_sequences:
                        fna_file.write(f">{query_id}__vs__{row_dict['subject_id']}\n")
                        fna_file.write(f"{fasta_sequences[row_dict['subject_id']]}\n")
                    else:
                        raise FileNotFoundError(f"Sequence {row_dict['subject_id']} not found in {query_id}_anc_ORFs.fna")
                    if row_dict['subject_id'] in translated_sequences:
                        faa_file.write(f">{query_id}__vs__{row_dict['subject_id']}\n")
                        faa_file.write(f"{translated_sequences[row_dict['subject_id']]}\n")
                    else:
                        raise FileNotFoundError(f"Translated sequence {row_dict['subject_id']} not found in {query_id}_anc_ORFs.faa")
                    written_pairs.add(pair)
        
        if mode == 'best' and best_hit:
            pair = (best_hit['query_id'], best_hit['subject_id'])
            if pair not in written_pairs:
                query_id = best_hit['query_id']
                fasta_sequences = read_fasta(f"{sanitize_filename(query_id)}_anc_ORFs.fna")
                translated_sequences = read_fasta(f"{sanitize_filename(query_id)}_anc_ORFs.faa")
                if best_hit['subject_id'] in fasta_sequences:
                    fna_file.write(f">{best_hit['query_id']}__vs__{best_hit['subject_id']}\n")
                    fna_file.write(f"{fasta_sequences[best_hit['subject_id']]}\n")
                else:
                    raise FileNotFoundError(f"Sequence {best_hit['subject_id']} not found in {query_id}_anc_ORFs.fna")
                if best_hit['subject_id'] in translated_sequences:
                    faa_file.write(f">{best_hit['query_id']}__vs__{best_hit['subject_id']}\n")
                    faa_file.write(f"{translated_sequences[best_hit['subject_id']]}\n")
                else:
                    raise FileNotFoundError(f"Translated sequence {best_hit['subject_id']} not found in {query_id}_anc_ORFs.faa")
                written_pairs.add(pair)

def get_ancorfs_fasta(tsv_files, mode):
    written_pairs = set()
    with open(f'ancORFs_{mode}.fna', 'w') as fna_file, open(f'ancORFs_{mode}.faa', 'w') as faa_file:
        for tsv in tsv_files:
            process_tsv_file(tsv, mode, fna_file, faa_file, written_pairs)

def main():
    parser = argparse.ArgumentParser(description='Build FASTA files with ancestral ORFs.')
    parser.add_argument("ssearch36_tsv", nargs='+', help="The ssearch36 output TSV files.")
    parser.add_argument('-m', '--mode', action='append', help='The filtering mode(s)')

    args = parser.parse_args()

    for mode in args.mode:
        if mode != 'best':
            try:
                float(mode)
            except ValueError:
                raise ValueError(f"Invalid mode: {mode}. Mode should be 'best' or a float value.")
        
        get_ancorfs_fasta(args.ssearch36_tsv, mode)

if __name__ == '__main__':
    main()