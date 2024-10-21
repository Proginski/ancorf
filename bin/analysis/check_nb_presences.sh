#!/bin/bash

grep -f <(sed "s/$/_elongated/" ../DENSE/table1_data/NCOS_CDS.txt) ../DENSE/dense/TRG_match_matrix.tsv | awk -F"\t" '{counter=0 ; for(i=3;i<=NF;i++){if($i != "noM"){counter++}} ; print counter}' | sort | uniq -c
