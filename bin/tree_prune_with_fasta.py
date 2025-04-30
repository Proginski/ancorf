#!/usr/bin/env python3

import argparse
from Bio import AlignIO
from Bio import Phylo

def prune_tree(input_tree_file, fasta_file, output_tree_file):
    # Read the input tree
    input_tree = Phylo.read(input_tree_file, "newick")

    # Read the FASTA file to get the list of taxa
    fasta = AlignIO.read(fasta_file, "fasta")
    fasta_taxa = [record.id for record in fasta]

    # Prune the tree to include only the taxa in the FASTA file
    taxa_to_remove = [taxon for taxon in input_tree.get_terminals() if taxon.name not in fasta_taxa]
    for taxon in taxa_to_remove:
        input_tree.prune(taxon)

    # Write the pruned tree to a file
    Phylo.write(input_tree, output_tree_file, "newick")

def main():
    parser = argparse.ArgumentParser(description="Prune a tree to match the taxa in a FASTA file.")
    parser.add_argument("input_tree", help="Input tree file in Newick format")
    parser.add_argument("fasta", help="FASTA file")
    parser.add_argument("output_tree", help="Output pruned tree file in Newick format")
    args = parser.parse_args()

    prune_tree(args.input_tree, args.fasta, args.output_tree)

if __name__ == "__main__":
    main()
