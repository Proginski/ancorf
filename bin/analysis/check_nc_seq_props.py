#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This program reads a newick tree.
For every _to_align.fna file in a directory, it reads the corresponding _to_align.tree file and write a fasta with the outgroup genomes.
"""

import argparse
import dendropy
from Bio import SeqIO
import glob
import sys
from Bio.SeqRecord import SeqRecord


def main():
    parser = argparse.ArgumentParser(description='Find MRCA of focal genome and closest outgroup genome.')
    parser.add_argument('tree_file', help='The Newick tree file of the genomes.')
    parser.add_argument('directory', help='The FASTA file with genome tags.')
    parser.add_argument('--first', '-f', action='store_true', help='Only use the first outgroup.')
    args = parser.parse_args()

    # print("Reading tree file...")
    tree = dendropy.Tree.get(path=args.tree_file,
        schema="newick",
        rooting='force-rooted')
    # print(tree.as_ascii_plot())

    # Get the list of file whose name ends with "_toalign_AA.fna" in the directory
    fasta_list = glob.glob(args.directory+'/*_toalign_AA.fna')
    # Take any file as an example
    fasta_file = fasta_list[0]
    focal = None
    for record in SeqIO.parse(fasta_file, "fasta"):
            tag = record.description.split()[1]
            if tag == "focal":
                focal = record.id
                print(f"Focal genome: {focal}", file=sys.stderr)
                break
    

    #### TREE INFORMATION ####
    # Open the tree file to order the neighbors with respect to the focal genome.
    # To do so, get the distance between the focal genome and each neighbor.
    # Two of more neighbors of the same "outgroup" should have the same distance to the focal genome.
    # The classical 'time of divergence' may lead to small differences if branches are not of equal length.
    # For example, Pan paniscus and Pan troglodytes form a possible "outgroup" in the human-and-its-neighbors tree,
    # but their calculated time of divergence with human may not be the same.
    # So just get the rank of the most recent common ancestor (MRCA) node.

    # Get the list of taxon labels (e.g. genome names)
    taxon_labels = [ tree.taxon_namespace[i].label for i in range(0,len(tree.taxon_namespace)) ]

    # Set every edge length to 1. By doing so, we are sure the calc_node_root_distances can be used
    # Besides, we do not care about the actual value of each edge length in this case.
    for node in tree.preorder_node_iter():
        node.edge.length = 1

    # Adds attribute “root_distance” to each node, with value set to the sum of 
    # edge lengths from the node to the root. Returns list of distances. 
    tree.calc_node_root_distances(return_leaf_distances_only=False) 

    # For each genome in taxon_labels, get its most recent common ancestor with 
    # the focal genome (= a node), and get its distance to the tree root.
    root_distance = { name:int(tree.mrca(taxon_labels=[focal, name]).root_distance) for name in taxon_labels if name != focal }
    root_distance[focal]=max(root_distance.values()) +1

    # Reverse the dictionary to make further comparisons more intuitive (focal is 0, its closest neighbor is 1, etc.)
    max_distance = max(root_distance.values())
    focal_distance = { genome: max_distance - root_distance[genome] for genome in root_distance }

    # List of genomes sorted by distance to the focal genome
    sorted_genomes = sorted(focal_distance, key=focal_distance.get)




    # For each FASTA :
    # print(f"Found {len(fasta_list)} files.")
    for fasta_file in fasta_list:
        genomes_names = []
        cds_names = []
        for record in SeqIO.parse(fasta_file, "fasta"):
            tag = record.description.split()[1]
            if tag == "genome":
                genomes_names.append(record.id)
            elif tag == "CDS":
                cds_names.append(record.id)

        # Get the maximum distance for a CDS match
        # Genome matches that are farther than this distance are "outgroup" matches
        CDS_distances = [focal_distance[name] for name in sorted_genomes if name in cds_names]
        genome_distances = [focal_distance[name] for name in sorted_genomes if name in genomes_names]

        if len(CDS_distances) > 0:
            max_CDS_distance = max(CDS_distances)
        else:
            max_CDS_distance = 0

        outgroup_genomes = [name for name in genomes_names if focal_distance[name] > max_CDS_distance]

        if args.first and len(outgroup_genomes) > 0:
            first_outgroup_distance = min([distance for distance in genome_distances if distance > max_CDS_distance])
            outgroup_genomes = [name for name in outgroup_genomes if focal_distance[name] == first_outgroup_distance]
        
        if len(outgroup_genomes) > 0:
            # Get the id of the sequences of the outgroup genomes
            outgroup_ids = [record.id for record in SeqIO.parse(fasta_file, "fasta") if record.id in outgroup_genomes]
            # Retrieve the sequence in the nucleotides alignment file
            nuc_aligned_file = fasta_file.replace("_toalign_AA.fna", "_aligned.fna")
            outgroup_records = [record for record in SeqIO.parse(nuc_aligned_file, "fasta") if record.id in outgroup_ids]
            # Convert this alignment fasta file into a regular fasta file (remove gaps in the sequences)
            for record in outgroup_records:
                record.seq = record.seq.ungap('-')
            # Write the fasta file 
            outgroup_file_nuc = fasta_file.replace("_toalign_AA.fna", "_outgroups.fna")
            SeqIO.write(outgroup_records, outgroup_file_nuc, "fasta")
            # SeqIO.write(outgroup_records, sys.stdout, "fasta")






if __name__ == '__main__':
    main()