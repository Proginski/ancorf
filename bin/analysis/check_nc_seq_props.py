#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
The goal of this program is to retrieve the non-coding sequences of the first outgroup genomes for each ancestral reconstruction.

This program reads a newick tree.
For every _to_align.fna file in a directory, it reads the corresponding _to_align.tree file and writes a fasta with the outgroup genomes.
"""

import argparse
import dendropy
from Bio import SeqIO
import glob
import os
import sys
from Bio.SeqRecord import SeqRecord


def main():
    parser = argparse.ArgumentParser(description='Find MRCA of focal genome and closest outgroup genome.')
    parser.add_argument('tree_file', help='The Newick tree file of the genomes.')
    parser.add_argument('input_directory', help='The input directory with the "_to_align" FASTA files with genome tags.')
    parser.add_argument('output_directory', help='The output directory to save the outgroup FASTA files.')
    parser.add_argument('--first', '-f', action='store_true', help='Only use the first outgroup.')
    args = parser.parse_args()

    # Ensure the output directory exists
    if not os.path.exists(args.output_directory):
        os.makedirs(args.output_directory)

    # Read the tree file
    tree = dendropy.Tree.get(path=args.tree_file, schema="newick", rooting='force-rooted')

    # Get the list of files whose name ends with "_toalign.fna" in the input directory
    fasta_list = glob.glob(os.path.join(args.input_directory, '*_toalign.fna'))

    # Take any file as an example in order to identify the focal genome
    fasta_file = fasta_list[0]
    focal = None
    for record in SeqIO.parse(fasta_file, "fasta"):
        tag = record.description.split()[1]
        if tag == "focal":
            focal = record.id
            print(f"Focal genome: {focal}", file=sys.stderr)
            break

    # Get the list of taxon labels (e.g., genome names)
    taxon_labels = [tree.taxon_namespace[i].label for i in range(len(tree.taxon_namespace))]

    # Set every edge length to 1
    for node in tree.preorder_node_iter():
        node.edge.length = 1

    # Calculate root distances
    tree.calc_node_root_distances(return_leaf_distances_only=False)

    # Get root distances for each genome
    root_distance = {name: int(tree.mrca(taxon_labels=[focal, name]).root_distance) for name in taxon_labels if name != focal}
    root_distance[focal] = max(root_distance.values()) + 1

    # Reverse the dictionary for further comparisons
    max_distance = max(root_distance.values())
    focal_distance = {genome: max_distance - root_distance[genome] for genome in root_distance}

    # List of genomes sorted by distance to the focal genome
    sorted_genomes = sorted(focal_distance, key=focal_distance.get)

    # Process each FASTA file
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
            outgroup_records = [record for record in SeqIO.parse(fasta_file, "fasta") if record.id in outgroup_ids]
            # Use the fasta name as the description
            for record in outgroup_records:
                record.description = fasta_file.split("/")[-1].replace("_toalign.fna", "")
            # Write the fasta file to the output directory
            if args.first:
                outgroup_file_nuc = os.path.join(args.output_directory, os.path.basename(fasta_file).replace("_toalign.fna", "_first_outgroup.fna"))
            else:
                outgroup_file_nuc = os.path.join(args.output_directory, os.path.basename(fasta_file).replace("_toalign.fna", "_outgroups.fna"))
            SeqIO.write(outgroup_records, outgroup_file_nuc, "fasta")


if __name__ == '__main__':
    main()