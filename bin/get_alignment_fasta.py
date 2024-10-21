#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This program generate a FASTA file with the CDS of the focal genome, and the CDS or the intergenic region of its neighbors.

It takes as input : 
- a text file (one column) with the list of CDS from the focal to align.
- a 'TRG_table.tsv' file from DENSE.
- a text file (one column) with the names of the genomes.
- the name of the focal genome.
"""


import argparse
from collections import defaultdict
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from pybedtools.bedtool import BedTool


parser = argparse.ArgumentParser(description='Generate a FASTA file with the CDS of the focal genome, and the CDS or the intergenic region of its neighbors.')
parser.add_argument('focal_list', type=str, help='Text file with the list of CDS from the focal to align.')
parser.add_argument('trg_table', type=str, help='TRG_table.tsv file from DENSE.')
parser.add_argument('names', type=str, help='txt file with the names of the genomes.')
parser.add_argument('focal', type=str, help='Name of the focal genome.')
args = parser.parse_args()

queries=[]
with open(args.focal_list) as fl:
    for line in fl:
        queries.append(line.rstrip())
names=[]
with open(args.names) as ns:
    for line in ns:
        names.append(line.rstrip())
focal=args.focal
print(f"Focal genome: {focal}")
print(f"Names of the genomes: {names}")

# Store the hits in a dictionary
hits_dic={}
neighbors_dic=defaultdict(lambda: defaultdict(list))

# For each line of the TRG table, if the query is in the list of queries, store the neighbors' hits in a dictionary
with open(args.trg_table) as tt:
    for line in tt:
        line=line.rstrip().split("\t")
        add_to_neighbors=True

        query=line[1].replace("_elongated","") # Remove the "_elongated" suffix from line[1]
        neighbor=line[2]
        
        target=line[3].replace("_elongated","")
        target_type=line[4]

        if query in queries and neighbor != focal: # exclude self-matches
            if query not in hits_dic:
                hits_dic[query] = {}
            if neighbor not in hits_dic[query]:
                hits_dic[query][neighbor] = {}
            if target_type not in hits_dic[query][neighbor]:
                if target_type == "genome" and "CDS" in hits_dic[query][neighbor]:
                    add_to_neighbors=False
                    continue
                else:
                    hits_dic[query][neighbor][target_type] = target
            else:
                sys.exit(f"Error: target_type already in the dictionary. {line}")
            if add_to_neighbors:
                neighbors_dic[neighbor][target_type].append(target)


# Build a nucleotide FASTA for the queries
entries_dic = defaultdict(lambda: defaultdict(list))

with open(focal+"_CDS_elongated.fna", "r") as input_handle:
    for record in SeqIO.parse(input_handle, "fasta"):
        record.id = record.id.replace("_elongated","")
        if record.id in queries:
            # Store the record in the dictionary
            entries_dic[focal][record.id] = record

# For each neighbor, build a nucleotide FASTA with the CDS collected in the dictionary
for name in names:
    if name != focal:
        print(name)
        # Build a nucleotide FASTA with the CDS targets
        with open(name+"_CDS_elongated.fna", "r") as input_handle: #, open(name+"_CDS_targets.fna", "w") as output_handle:
            for record in SeqIO.parse(input_handle, "fasta"):
                record.id = record.id.replace("_elongated","")
                if record.id in neighbors_dic[name]["CDS"]:
                    entries_dic[name][record.id] = record

# Build a nucleotide FASTA for further alignment
# Only keep queries with at least --num_outgroups genome (i.e. non-coding) outgroups neighbors 
for query in hits_dic:
    if len(hits_dic[query]) > 1:
        print(query)
        # Build a nucleotide FASTA file with the CDS of the focal genome, and the CDS or the intergenic region of its neighbors
        def sanitize_filename(filename):
            return filename.replace(":", "__COLON__").replace("/", "__SLASH__").replace("!", "__EXCLAMATION__").replace("|", "__PIPE__")

        query_sanitized = sanitize_filename(query)
        with open(f"{query_sanitized}_toalign.fna", "w") as output_handle:
            # Change the header to the name of the genome
            print(entries_dic[focal][query])
            entries_dic[focal][query].id = focal
            entries_dic[focal][query].description = "focal"
            SeqIO.write(entries_dic[focal][query], output_handle, "fasta")
            # Also write the CDS of the neighbors
            for neighbor in hits_dic[query]:
                if "CDS" in hits_dic[query][neighbor]:
                    # Change the header to the name of the genome
                    entries_dic[neighbor][hits_dic[query][neighbor]["CDS"]].id = neighbor
                    entries_dic[neighbor][hits_dic[query][neighbor]["CDS"]].description = "CDS"
                    SeqIO.write(entries_dic[neighbor][hits_dic[query][neighbor]["CDS"]], output_handle, "fasta")
                else: # genome type
                    target=hits_dic[query][neighbor]["genome"]
                    target=target.split("_")
                    seq="_".join(target[:-2])
                    start=int(target[-2])
                    end=int(target[-1])
                    strand = "+"
                    if start > end:
                        strand = "-"
                        start,end = end,start
                    start -= 1
                    # Make a bedtool object with the intergenic region
                    bed = BedTool(f"{seq}\t{start}\t{end}\t.\t.\t{strand}", from_string=True)
                    # Then, use bedtools getfasta to extract its sequence
                    fna = bed.sequence(fi=neighbor+".fna",s=True)
                    # Change the header to the name of the genome
                    for record in SeqIO.parse(fna.seqfn, "fasta"):
                        record.id = neighbor
                        record.description = "genome"
                        SeqIO.write(record, output_handle, "fasta")