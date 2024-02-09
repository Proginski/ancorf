#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from ete3 import Tree

def prune_tree(newick_file, taxa_file, output_file=None):
    with open(taxa_file, 'r') as f:
        taxa_to_keep = [line.strip() for line in f]

    tree = Tree(newick_file)
    tree.prune(taxa_to_keep)

    if output_file:
        tree.write(outfile=output_file, format=1)
    else:
        print(tree.write(format=1))

def main():
    parser = argparse.ArgumentParser(description='Prune a Newick tree.')
    parser.add_argument('newick_file', help='The Newick file to prune.')
    parser.add_argument('taxa_file', help='The file with the list of taxa to keep.')
    parser.add_argument('-o', '--output', help='The output file name.')

    args = parser.parse_args()

    prune_tree(args.newick_file, args.taxa_file, args.output)

if __name__ == '__main__':
    main()