#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The aim of the script is to detect abberrant trees when it comes to reconstruct the last common ancestor between a focal leaf and the first non-coding outgroup. 
This script takes as input a focal leaf, and a TSV file with orf\tnewick_string\tancestor_node\tCDS_leaves.
It writes a new TSV file with orf\tboolean\tancestor_node\toutgroup_CDS_leaves.
"""

import argparse
from ete3 import Tree

def get_outgroup_leaves(tree, node_name, focal_leaf_name):
    node = tree.search_nodes(name=node_name)
    if not node:
        raise ValueError(f"Node with name '{node_name}' not found in the tree.")
    node = node[0]  # Assuming the node name is unique

    focal_leaf = tree.search_nodes(name=focal_leaf_name)
    if not focal_leaf:
        raise ValueError(f"Focal leaf with name '{focal_leaf_name}' not found in the tree.")
    focal_leaf = focal_leaf[0]  # Assuming the focal leaf name is unique

    # Traverse from the focal leaf up to the input node
    current_node = focal_leaf
    path = []
    while current_node and current_node != node:
        path.append(current_node.name)
        current_node = current_node.up

    if current_node == node:
        path.append(node.name)
    else:
        raise ValueError(f"Focal leaf '{focal_leaf_name}' is not a descendant of node '{node_name}'.")

    # CDS_node = path[-2]

    # # Get all leaves under the input node
    # all_leaves = set(node.get_leaf_names())

    # CDS_leaves = set(tree.search_nodes(name=CDS_node)[0].get_leaf_names())
  
    # outgroup_leaves = all_leaves - CDS_leaves

    # return outgroup_leaves

    # The provided "node_name" has two children branches (bifurcation).
    # One contains the focal leaf (plus possibly other nodes and leaves), and the other does not.
    # The "penultimate_node" is the node that contains all the leaves of this first branch (the focal's branch)..
    penultimate_node = path[-2]

    # All the leaves of the tree
    all_leaves = set(tree.get_leaf_names())

    # All the leaves of the penultimate node
    penultimate_node_leaves = set(tree.search_nodes(name=penultimate_node)[0].get_leaf_names())

    other_leaves = all_leaves - penultimate_node_leaves

    return other_leaves

def main():
    parser = argparse.ArgumentParser(description="Print all nodes from the focal leaf to the input node in a Newick tree.")
    parser.add_argument("input_tsv", help="Path to the input TSV file.")
    parser.add_argument("output_tsv", help="Path to the output TSV file.")
    parser.add_argument("-f", "--focal_leaf_name", help="Name of the focal leaf.")
    args = parser.parse_args()

    with open(args.input_tsv, 'r') as infile, open(args.output_tsv, 'w') as outfile:
        infile.readline()  # Read and skip the header line of the input file
        header = "orf\ttreeIsOK\tanode\tpb_leaves\n"
        outfile.write(header)
        for line in infile:
            orf, newick_string, ancestor_node, cds_leaves = line.strip().split('\t')
            if cds_leaves == '.':
                outfile.write(f"{orf}\tTRUE\t{ancestor_node}\t.\n")
            else:
                try:
                    tree = Tree(newick_string, format=1)
                    outgroup_leaves = get_outgroup_leaves(tree, ancestor_node, args.focal_leaf_name)
                    # The problematic leaves are those which are outgroup and CDS
                    outgroup_CDS_leaves = set(outgroup_leaves) & set(cds_leaves.split('|'))
                    outfile.write(f"{orf}\t{str(not bool(outgroup_CDS_leaves)).upper()}\t{ancestor_node}\t{'|'.join(outgroup_CDS_leaves)}\n")
                except ValueError as e:
                    outfile.write(f"{orf}\tFALSE\t{ancestor_node}\tERROR: {str(e)}\n")
                except Exception as e:
                    outfile.write(f"{orf}\tFALSE\t{ancestor_node}\tERROR: {str(e)}\n")

if __name__ == "__main__":
    main()