#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 17 16:40:15 2022

@author: paul.roginski
"""


import argparse
import dendropy
import pandas as pd
import sys


# Arguments parsing
parser = argparse.ArgumentParser()
parser.add_argument("-names", required=True, help="csv file matching names (col1) and phylip_names (col2)", default=False)
parser.add_argument("-tree", required=True, help="newick file for the input phylogeny tree")
parser.add_argument("-out", required=True, help="output newick file")
args = parser.parse_args()
    
# newick file for the phylogeny tree
tree_file = args.tree

names = pd.read_csv(args.names, names = ["name", "phylip_name"], sep = "\t")

tree = dendropy.Tree.get(path=tree_file, schema='newick', preserve_underscores=True)

# Only keep taxons that are in the names file
for taxon in tree.taxon_namespace :
    if taxon.label not in names["name"].tolist() :
        tree.taxon_namespace.remove_taxon(taxon)

new_names = []
for name in names["name"] :
    node_to_change = tree.find_node_with_taxon_label(name)
    new_name = names.loc[names['name'] == name, 'phylip_name'].item()
    
    if new_name in new_names or (new_name != name and new_name in names["name"]):
        sys.exit("Error : the new phylip name " + new_name + " is already in the tree. The genome names must be changed to differ even when shortened to phylip format.")
    else :
        node_to_change.taxon.label = new_name

tree.write(path=args.out, schema="newick")