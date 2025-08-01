# Ancestral Open Reading Frame Reconstruction Pipeline for *de novo* genes

A Nextflow workflow for identifying ancestral ORFs that gave birth to present-day de novo genes through phylogenetic reconstruction and sequence analysis.

## Overview

This NEXTFLOW workflow takes a list of de novo emerged genes in a focal genome, along with their homologous sequences/regions across neighbor genomes. Using PRANK or PREQUEL, it reconstructs an ancestral DNA locus that is assumed to be non-coding. From this locus it extracts all the ORFs (stop-to-stop) >= 20 codons, and aligns them to the current de novo CDS. It thus identifies ancestral ORFs that gave birth to today's de novo genes.

## Input Parameters

- **focal** = the name of the focal genome
- **gendir** = the directory with the genomes GFF3 and FASTA files (focal and neighbors)
- **tree** = a NEWICK tree with the focal and neighbor genomes
- **TRG_table** = a DENSE workflow output (TSV file) with the precomputed matches of the de novo CDS with CDS or genome regions of the neighbor genomes
- **queries** = a text file with the desired list of de novo CDSs to process
- **mode** = different modes of ancestral locus reconstruction (default: 'prank')
- **outdir** = name of the results directory

## Usage

```
nextflow run proginski/ancorf -profile <SINGULARITY-APPTAINER-DOCKER> --gendir <DIR_WITH_GFF_AND_FASTA> --focal <FOCAL_GENOME_NAME> --tree <NEWICK_WITH_FOCAL_AND_NEIGHBORS> --TRG_table <TRG_TABLE> --queries <QUERIES> --mode <PRANK-PREQUEL> --outdir <OUTDIR>
```

## Container Requirements

The workflow is expected to be run with a container manager (Singularity, Apptainer, Docker). It automatically pulls the DockerHub image proginski/ancorf.

## Workflow Architecture

- **main.nf** is just a generic entry to the specific workflows/ancorf.nf workflow
- **workflows/ancorf.nf** orchestrates the execution of the different modules
- **modules/local/ancorf_modules.nf** contains the different modules
- **bin/** contains the utility scripts
- **containers/** contains the pulled images
- **nextflow.config** the general config file
- **nextflow_schema.json** the helper json for the parameters

## Processes

### Core Modules

- **CHECK_INPUTS**: check the integrity of the gendir, tree, and TRG_table
- **EXTRACT_CDS**: extract the CDSs from the genome GFF3 and FASTA files
- **ELONGATE_CDS**: elongate the CDSs of 99 nucleotides upstream and downstream to add local genomic context for the reconstruction
- **ALIGNMENT_FASTA**: build an alignment FASTA file for each query de novo CDS with its CDS or genome matches across the neighbor genomes

### Ancestral ORF Reconstruction

#### ANCESTRAL_ORFS_PRANK (mode: "prank")
Align the entries of the alignment FASTA file with macse_v2.07, fix the format and run prank providing the genomes' tree without iteration (-once), preserving the input alignment (-keep), and providing the different inferred ancestors in the tree (-showanc). It then selects the most recent common ancestor between the last CDS match, and the first outgroup (see https://github.com/i2bc/dense). It extracts all the ORFs (stop to stop) of at least 20 codons across the three frames. It finally tries to align the query de novo CDS with these ancestral ORFs, using ssearch36 (e-value 1e-2).

#### ANCESTRAL_ORFS_PREQUEL (mode: "prequel")
Same as ANCESTRAL_ORFS_PRANK, except that for the reconstruction of the ancestral locus, it uses PREQUEL. It first builds a model with phyloFit based on the alignment and the provided genomes tree. It preserves the input tree in the model that is provided to prequel.

### Output Generation

#### ANCORFS_FASTA
Build FASTA files with the selected ancestral ORFs:
- **best**: only the best matching ORF (e-value based) for each query
- **1e-2**: all matches < 1e-2
- **1e-3**: all matches < 1e-3

## Dependencies

DENSE workflow: https://github.com/i2bc/dense

## Citation

If you use this pipeline in your research, please cite the associated publication and the DENSE workflow.
