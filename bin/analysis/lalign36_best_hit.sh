#!/bin/bash

# ancorf_tree_genomes/rna-gnl__PIPE__WGS__COLON__LPNL__PIPE__AWRI3578_g3804_mrna_anc_ORFs_vs_rna-gnl__PIPE__WGS__COLON__LPNL__PIPE__AWRI3578_g3804_mrna_lalign36.tsv
awk 'BEGIN{FS=OFS="\t"} !/^#/ {print gensub(/_vs_.*/,"","g",FILENAME),$2; exit}' $1 | while read rawname_ORF
do
	rawname=$(echo "${rawname_ORF}" | awk -F"\t" '{print $1}')
	ORF=$(echo "${rawname_ORF}" | awk -F"\t" '{print $2}')

	faOneRecord ${rawname}.faa $ORF
done
