process GET_CANDIDATES_LIST {

	publishDir 'ANC_ORF'
	
	input:
		path candidates

	output:
		path "candidates.txt"

	"""	
	#cat $candidates | shuf -n 30 > candidates.txt
	cp $candidates candidates.txt
	
	# Since the query's name will be used to name most of the subsequent files,
	# some special characters in it can stop the pipeline before the end.
	# Correct it now.
	sed -i "s/:/__COLON__/g" candidates.txt

	if [ ! -s candidates.txt ]; then echo "candidates.txt DOES NOT EXISTS OR IS EMPTY."; exit 1; fi
	"""
}


process CHECK_INPUTS {

	publishDir 'ANC_ORF'
	
	input:
		path original_tree
		path species_list
		val focal
		path candidates

	output:
		path "neighbors.txt", emit : neighbors
		path "${focal}_CDS_multielongated_candidates.faa", emit : multicandidates_faa
		
	"""
	# Check the existence of the following files :
	for file in $original_tree $species_list; do
		if [ -s \$file ]; then
			echo "\$file exists and is not empty."
		else
			echo "\$file DOES NOT EXISTS or is EMPTY."
			exit 1
		fi
	done
	
	
	# Check that the focal species is indeed in the species list.
	if grep -q -x $focal $species_list; then
		echo "$focal is found in $species_list"
	else
		echo "$focal CANNOT BE FOUND in $species_list"
		exit 1	
	fi
	
	
	# Get the species list with only the neighbors (remove the focal)
	grep -v -x $focal $species_list > neighbors.txt
	
	
	# If necessary, generate the FASTA with the CDS of the focal species, elongate by 100 nucleotides up and downstream and such translated :
	
	# 100_nucl_translated_in_frame_0____translated_CDS____100_nucl_translated_in_frame_0
	# 100_nucl_translated_in_frame_0____translated_CDS____100_nucl_translated_in_frame_1
	# 100_nucl_translated_in_frame_0____translated_CDS____100_nucl_translated_in_frame_2
	# 100_nucl_translated_in_frame_1____translated_CDS____100_nucl_translated_in_frame_0
	# 100_nucl_translated_in_frame_1____translated_CDS____100_nucl_translated_in_frame_1
	# 100_nucl_translated_in_frame_1____translated_CDS____100_nucl_translated_in_frame_2
	# 100_nucl_translated_in_frame_2____translated_CDS____100_nucl_translated_in_frame_0
	# 100_nucl_translated_in_frame_2____translated_CDS____100_nucl_translated_in_frame_1
	# 100_nucl_translated_in_frame_2____translated_CDS____100_nucl_translated_in_frame_2
	
	# This will be further refered as the multielongated translated CDS FASTA.
	
	if [ ! -s ${PWD}/GENOMES/${focal}_CDS_multielongated.faa ]
	then /SCRIPTS/get_CDS_borders.sh --gfasta ${PWD}/GENOMES/${focal}.fna --gff ${PWD}/GENOMES/${focal}.gff --genome ${PWD}/GENOMES/${focal}.genome --size 100 --multiframe True --outdir ${PWD}/GENOMES/
	fi
	
	cat $candidates | sed "s/^/>/" | sed "s/\$/_elongated/" | grep -f - ${PWD}/GENOMES/${focal}_CDS_multielongated.faa -A1 > ${focal}_CDS_multielongated_candidates.faa
	"""
}


process CHECK_NEIGHBORS {

	clusterOptions  "-q lowprio -l ncpus=40 -l walltime=43800:00:00"

	input:
		val neighbor
		val focal
		path multicandidates_faa
	
	output:
		val true
		
	"""
	# In order to search for homologous CDS in the neighbor, get a FASTA with the elongated CDS of the neighbor (100 nucl upstream and downstream).
	# In case of multi matchs, the one with the best evalue is more likely to be the right one with this elongated version (bring genomic context).
	if [ ! -s ${PWD}/GENOMES/BLAST_DB/${neighbor}_CDS_elongated.faa.ptf ]
	then
		# Generate the elongated translated CDS FASTA.
		bash /SCRIPTS/get_CDS_borders.sh --gfasta ${PWD}/GENOMES/${neighbor}.fna --gff ${PWD}/GENOMES/${neighbor}.gff --genome ${PWD}/GENOMES/${neighbor}.genome --size 100 --multiframe False --outdir ${PWD}/GENOMES/
		# Create the blast db.
		makeblastdb -dbtype prot -in ${PWD}/GENOMES/${neighbor}_CDS_elongated.faa -out ${PWD}/GENOMES/BLAST_DB/${neighbor}_CDS_elongated.faa
	fi
	
	
	# Also create a nucleotide db with the entier genome, to look for any matchs (only considered if no CDS matchs are found).
	if [ ! -s ${PWD}/GENOMES/BLAST_DB/${neighbor}.fna.ntf ]
	then
		makeblastdb -dbtype nucl -in ${PWD}/GENOMES/${neighbor}.fna -out ${PWD}/GENOMES/BLAST_DB/${neighbor}.fna
	fi 
	
	
	# If necessary create the BLASTP directory.
	mkdir -p ${PWD}/BLASTP/
	
	
	# If not already done, perform a BLASTp of the multielongated translated CDS FASTA of the focal species against the elongated CDS FASTA of its neighbor.
	# The -num_threads option should be adapted.
	if [ ! -s ${PWD}/BLASTP/${focal}_CDS_multielongated_candidates_blastp_${neighbor}_CDS_elongated.out ]
	then
		blastp -query $multicandidates_faa -db ${PWD}/GENOMES/BLAST_DB/${neighbor}_CDS_elongated.faa -num_threads 40 -outfmt "7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen" -evalue 0.001 > ${PWD}/BLASTP/${focal}_CDS_multielongated_candidates_blastp_${neighbor}_CDS_elongated.out
	fi 
	
	
	# Similarly, perform if nedded, a tBLASTn of the multielongated translated CDS FASTA of the focal species against the genomic nucleotide FASTA of its neighbor.
	if [ ! -s ${PWD}/TBLASTN/${focal}_CDS_multielongated_candidates_tblastn_${neighbor}.out ]
	then
		tblastn -query $multicandidates_faa -db ${PWD}/GENOMES/BLAST_DB/${neighbor}.fna -num_threads 40 -outfmt "7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen" -evalue 0.001 > ${PWD}/TBLASTN/${focal}_CDS_multielongated_candidates_tblastn_${neighbor}.out
	fi 
	"""
}


process GET_PHYLIP_NAMES {

	publishDir 'ANC_ORF'
	
	input:
		path species_list
		path original_tree

	output:
		path "${species_list}_phylip_names", emit : phylip_names
		path "phylip_names_tree.nwk", emit : phylip_names_tree

	"""
	# Get a two columns csv file with original species names as first col and PHYLIP compatible names as the second col.
	awk '{short=substr(\$0,1,5)NR; lim_short=substr(short,1,10); final=gensub(/_/,"","g",lim_short); print \$0","final }' $species_list > ${species_list}_phylip_names
	
	# If duplicate taxon label exist in the input newick tree, if they are not in the species_list, correct them, otherwise prompt an error message
	bash /SCRIPTS/correct_duplicate_taxon_labels.sh --species_list $species_list --main_dir $PWD --process_dir \$PWD --tree $original_tree
	
	# Get a new tree with PHYLIP compatible names.
	python /SCRIPTS/get_phylip_named_tree.py -names ${species_list}_phylip_names -tree ${original_tree}_checked -out phylip_names_tree.nwk
	"""
}


process GET_ALIGNMENT_FASTA {

	clusterOptions  "-q bim -l ncpus=1 -l walltime=43800:00:00"

	publishDir 'ANC_ORF/FASTA_to_be_aligned'

	input:
		val check_neighbors_done
		val query
		val focal_species
		path species_list
		path phylip_names
		path tree

	output:
		// These output are optional since CDS with less than 2 neighbors with an homolog are filtered out at thise step.
		path "${query}_toalign.fna", emit : toali, optional: true
		path "${query}.nwk", emit : newick, optional: true
		path "${query}_names_n_types.txt", emit : names_n_types, optional: true
		
	"""
	bash /SCRIPTS/get_alignment_fasta.sh --query $query --species_list $species_list --main_dir $PWD --process_dir \$PWD --phylip_names $phylip_names --focal_species $focal_species --tree $tree
	"""
}


process NUCLEOTIDE_ALIGNMENT {

	clusterOptions  "-q bim -l ncpus=1 -l walltime=43800:00:00"

	publishDir 'ANC_ORF/FASTA_to_be_aligned/macse_results'

	input:
		path toali
		path newick
		path names_n_types

	output:
		path "*_aligned.fna", emit : NT_alignment
		path "*_aligned.faa", emit : AA_alignment
		path newick, emit : newick
		path names_n_types, emit : names_n_types
		
	"""
	# Align
	java -jar /packages/macse_v2.06.jar -prog alignSequences -seq $toali
	
	out_NT=\$(echo $toali | sed "s/_toalign.fna/_aligned.fna/")
	mv ${toali.baseName}_NT.fna \$out_NT
	out_AA=\$(echo $toali | sed "s/_toalign.fna/_aligned.faa/")
	mv ${toali.baseName}_AA.fna \$out_AA
	
	# Correct the NT file by replacing '!' characters by '-'
	sed -i 's|!|-|g' \$out_NT \$out_AA
	"""
}


process CONVERT_ALIGNED_FASTA_TO_PHYLIP {

	clusterOptions  "-q bim -l ncpus=1 -l walltime=43800:00:00"

	publishDir 'ANC_ORF/FASTA_to_be_aligned/PHYLIP'

	input:
		path NT_alignment
		path newick
		path names_n_types

	output:
		path "${NT_alignment.baseName}.phylip", emit : phylip
		path newick, emit : newick
		path NT_alignment, emit : NT_alignment
		path names_n_types, emit : names_n_types
		
	"""
	seqret -sequence FASTA::${NT_alignment} -outseq PHYLIP::${NT_alignment.baseName}.phylip
	"""
}


process PHYML_TREE {

	clusterOptions  "-q bim -l ncpus=1 -l walltime=43800:00:00"

	publishDir 'ANC_ORF/FASTA_to_be_aligned/PHYLIP'

	input:
		path phylip
		path newick
		path NT_alignment
		path names_n_types
		
	output:
		path "${phylip}_phyml_tree.txt", emit : phyml_tree
		path NT_alignment, emit : NT_alignment
		path names_n_types, emit : names_n_types
		
	"""
	### WARNING phyml needs a relative path after -i ###
	### WARNING depending on the computer used to install phyml 
	### (e.g. using nextflow with a singularity image built with another machine), 
	### launching phyml on another machine (or node, for HPCs) might prompt an "Illegal instruction" error, 
	### because CPUs are not configured the same way. 
	### See https://wiki.parabola.nu/Fixing_illegal_instruction_issues ###
	phyml -i $phylip -d nt -v e -o lr -c 4 -a e -b 0 -f e -u $newick
	"""
}


process PRANK {

	clusterOptions  "-q bim -l ncpus=1 -l walltime=43800:00:00"

	publishDir 'ANC_ORF/FASTA_to_be_aligned/PRANK'

	input:
		path NT_alignment
		path phyml_tree
		path names_n_types

	output:
		path "${NT_alignment.baseName}.best.anc.fas", emit : prank_fas
		path "${NT_alignment.baseName}.best.anc.dnd", emit : prank_tree
		path names_n_types, emit : names_n_types
		
	"""
	prank -d=${NT_alignment} -t=${phyml_tree} -o=${NT_alignment.baseName} -once -showanc
	"""
}


process GET_ANCESTOR_FNA {

	clusterOptions  "-q bim -l ncpus=1 -l walltime=43800:00:00"

	publishDir 'ANC_ORF/ANCESTORS'

	input:
		path prank_fas
		path prank_tree
		path phylip_names_tree
		path names_n_types

	output:
		path "*_anc.fna", emit : ancestor_fna, optional: true
		
	"""
	focal=\$(awk 'FNR == 1 {print \$1}' $names_n_types)
	IGR=\$(awk '\$2 == "IGR" {print \$1}' $names_n_types)
	name=\$(echo $names_n_types | sed "s/_names_n_types.txt\$//")
	
	# TEMP CORRECTION !
	if [ ! -z \$(echo "\${IGR:0:1}") ]; then
	
		python /SCRIPTS/tree_find_closest.py -tree $phylip_names_tree -focal \$focal -names <(echo "\${IGR}") -out \${name}_outgroup_names.txt
		
		python /SCRIPTS/tree_mrca.py -tree $prank_tree -names \${name}_outgroup_names.txt > \${name}_one_to_choose.txt
		
		echo ">\$(cat \${name}_one_to_choose.txt)" | grep -x -f - $prank_fas -A1 | sed "s/-//g" > \${name}_anc.fna
		
	fi
	"""
}


process GET_ANCESTOR_ORFS {

	clusterOptions  "-q bim -l ncpus=1 -l walltime=43800:00:00"

	publishDir 'ANC_ORF/ANCESTORS'

	input:
		path ancestor_fna
	output:
		path "${ancestor_fna.baseName}_ORFs.faa", emit : ancestor_ORF
	"""
	# -s 2 = stop to stop; -ml 60 = minimal length 60 nucl; -strand plus = positive strand only
	ORFfinder -in $ancestor_fna -s 2 -ml 60 -n false -strand plus > ${ancestor_fna.baseName}_ORFs.faa
	"""
}


process LALIGN {

	clusterOptions  "-q bim -l ncpus=1 -l walltime=43800:00:00"

	publishDir 'ANC_ORF/ANCESTORS'

	input:
		path ancestor_ORF
		val focal_species
	output:
		path "*.frags"
		path "*.faa"
	"""
	name=\$(echo $ancestor_ORF | sed "s/_anc_ORFs.faa\$//")
	# Retrieve the original name.
	name_uncorr=\$(echo \$name | sed "s/__COLON__/:/g")
	
	grep -x ">\${name_uncorr}" -A1 ${PWD}/GENOMES/${focal_species}_CDS.faa > \${name}.faa
	
	
	# Perform lalign without a tabular blast output format and and minimal e-value of 0.01
	# For each subject fragment, only keep the best (first) one.
	lalign36 -m8C -E 0.01 \${name}.faa $ancestor_ORF | \
	awk -F"\t" '\$0 !~ /^#/ && !(\$2 in data) {data[\$2] = 1; print \$0}' > \${name}.frags
	"""
}
