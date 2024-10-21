#!/bin/bash

# $1 input file
# $2 evalue
if [[ $2 == "" ]]; then	echo "Please provide an evalue treshold as second argument" ; fi

awk -F"\t" -v evalue="${2}" '!/^#/ && (!($2 in data)) && $11 < evalue {print $2 ; data[$2]=1}' $1
