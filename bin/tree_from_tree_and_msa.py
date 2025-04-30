#!/usr/bin/env python3

import argparse
from Bio import AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.BaseTree import Tree
import numpy as np
from scipy.optimize import minimize

def read_alignment(alignment_file):
    alignment = AlignIO.read(alignment_file, "fasta")
    return alignment

def read_tree(tree_file):
    tree = Phylo.read(tree_file, "newick")
    return tree

def calculate_pairwise_distances(alignment):
    calculator = DistanceCalculator('identity')
    dm = calculator.get_distance(alignment)
    return dm

def calculate_likelihood(tree, alignment, pairwise_distances):
    likelihood = 0.0
    for clade in tree.get_terminals():
        for other_clade in tree.get_terminals():
            if clade != other_clade:
                observed_distance = pairwise_distances[clade.name, other_clade.name]
                expected_distance = clade.branch_length + other_clade.branch_length
                likelihood -= (observed_distance - expected_distance) ** 2
    return likelihood

def optimize_branch_lengths(tree, alignment):
    pairwise_distances = calculate_pairwise_distances(alignment)

    def objective(branch_lengths):
        for i, clade in enumerate(tree.get_terminals()):
            clade.branch_length = branch_lengths[i]
        return -calculate_likelihood(tree, alignment, pairwise_distances)

    # Initialize branch lengths if they are None
    initial_branch_lengths = [clade.branch_length if clade.branch_length is not None else 0.1 for clade in tree.get_terminals()]
    result = minimize(objective, initial_branch_lengths, method='L-BFGS-B')
    optimized_branch_lengths = result.x

    for i, clade in enumerate(tree.get_terminals()):
        clade.branch_length = optimized_branch_lengths[i]

    return tree

def main():
    parser = argparse.ArgumentParser(description="Refine branch lengths of a tree based on an MSA without changing its topology.")
    parser.add_argument("alignment", help="Input alignment file in FASTA format")
    parser.add_argument("tree", help="Input tree file in Newick format")
    parser.add_argument("output_tree", help="Output refined tree file in Newick format")
    args = parser.parse_args()

    alignment = read_alignment(args.alignment)
    tree = read_tree(args.tree)
    refined_tree = optimize_branch_lengths(tree, alignment)
    Phylo.write(refined_tree, args.output_tree, "newick")

if __name__ == "__main__":
    main()