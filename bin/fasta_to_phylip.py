#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from Bio import AlignIO

# Create the parser
parser = argparse.ArgumentParser(description='Convert FASTA to PHYLIP format')

# Add the arguments
parser.add_argument('InputFile', metavar='input', type=str, help='the input FASTA file')
parser.add_argument('OutputFile', metavar='output', type=str, help='the output PHYLIP file')

# Parse the arguments
args = parser.parse_args()

# Convert the file
alignments = AlignIO.parse(args.InputFile, "fasta")
AlignIO.write(alignments, args.OutputFile, "phylip")

print(f"Converted {args.InputFile} to {args.OutputFile} in PHYLIP format.")