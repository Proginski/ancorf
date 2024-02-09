#!/bin/bash

A=$1
B=$2

# Write a header
awk 'BEGIN {OFS="\t" ; print "query","subject","len","score"}'

# In a loop, perform the pairwise alignments
awk '/>/ {print substr($1,2)}' ${A} | sed "s/>//" | while read ORF
do
	faOneRecord ${A} ${ORF} > A.faa

	water -asequence A.faa -bsequence ${B} water -gapopen 10.0 -gapextend 0.5 -outfile tmp -aformat score

	awk -v ORF="${ORF}" '
		BEGIN {OFS="\t"}
		FNR==1 {
			$1=ORF
			# Remove the parentheses from the score field
			gsub(/\(|\)/,"",$4)
			print $0
		}
	' tmp
done
