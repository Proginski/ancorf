#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script takes as input a tree file in newick format and a pair of leaf names.
It warns if the two leaves are not #!/usr/bin/env python3
"""

import argparse
import dendropy

# Parse command line arguments
parser = argparse.ArgumentParser(description="Check if two leaves are siblings in a tree.")
parser.add_argument("newick_files", nargs='+', help="The Newick files containing the trees.")
parser.add_argument("--leaf1", help="The name of the first leaf.")
parser.add_argument("--leaf2", help="The name of the second leaf.")
args = parser.parse_args()

print("Newick file\tareSiblings")

for newick_file in args.newick_files:
    try:
        tree = dendropy.Tree.get(path=newick_file, schema="newick")
        # rest of your code to process the tree
    except dendropy.DataParseError:
        print(f"Error: {newick_file} is not a valid Newick file.")
        continue

    # Find the nodes corresponding to the leaves
    node1 = tree.find_node_with_taxon_label(args.leaf1)
    node2 = tree.find_node_with_taxon_label(args.leaf2)

    # Check if the nodes were found in the tree
    proceed = True
    if node1 is None:
        print(f"Error: {args.leaf1} was not found in the tree. Available leaf names are: {', '.join([taxon.label for taxon in tree.taxon_namespace])}")
        proceed = False
    if node2 is None:
        print(f"Error: {args.leaf2} was not found in the tree. Available leaf names are: {', '.join([taxon.label for taxon in tree.taxon_namespace])}")
        proceed = False

    # Check if the nodes are siblings
    if proceed:
        if node1.parent_node != node2.parent_node:
            print(f"{newick_file}\tFALSE")
        else:
            print(f"{newick_file}\tTRUE")