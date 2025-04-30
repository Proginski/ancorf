process CHECK_INPUTS {

	input:
		path gendir
		path tree
		path TRG_table
		
	output:
		path "genome_files.txt",  emit : genome_files
		path "phylip_names.tsv", emit : phylip_names
		path "phylip_names_tree.nwk", emit : tree
		path TRG_table, emit : trg_table

	"""
	chmod -R +x ${projectDir}/bin

	# 1. Check that there is a single non-empty valid FASTA and GFF3 file per genome :

	# The list of non-empty files from 'gendir' :
	gendirlist=\$(find ${gendir}/. -maxdepth 1 -not -empty -ls | awk '{print \$NF}')

	echo "\${gendirlist}" | grep "\\.fas\$\\|\\.fna\$\\|\\.fasta\$\\|\\.fa\$" | sed -E "s~.*/(.*)\\..*~\\1~" | sort | uniq -c > nb_fasta.txt
	echo "\${gendirlist}" | grep "\\.gff\$\\|\\.gff3\$" | sed -E "s~.*/(.*)\\..*~\\1~" | sort | uniq -c > nb_gff.txt
	
	awk '
		# Record every genome 
		{ gn[\$2] = 1 }

		# For any line, if the nb of occurences is not one, record the genome in "pb"
		\$1 != 1 { pb[\$2]=1 }

		# Fasta list
		FNR==NR { fas[\$2] = 1 }
		# GFF list
		FNR!=NR { gff[\$2] = 1 }

		END {
			for(file in fas){ if (!(file in gff)) { pb[file] = 1 } }
			for(file in gff){ if (!(file in fas)) { pb[file] = 1 } }

			for(file in gn){ 
				if (file in pb) { print file >> "problematic_genomes.txt" }
				else { print file >> "valid_genomes.txt" } 
			}
		}
		' nb_fasta.txt nb_gff.txt
	
	if [ -s problematic_genomes.txt ]
	then
		echo "WARNING : Each genome must be provided with one genomic FASTA file ('.fna','.fa','.fas','.fasta') and one GFF3 annotation file ('.gff3','.gff')"
		echo "ERROR : some genomes seem to not have a single non-empty valid FASTA and GFF3 file (see 'problematic_genomes.txt', 'nb_fasta.txt', 'nb_gff.txt') :"
		cat problematic_genomes.txt
		exit 1
	fi

	# The list of the genomes to process is generated from the list of FASTA files in gendir.
	echo "Genomes to procces (based on FASTA files in $gendir):"
	cat valid_genomes.txt
	echo ""


	# 2. Check that GFF files are GFF3 compliant and have "CDS", "gene", and "mRNA" features.
	cat valid_genomes.txt | while read genome
	do
		gff=$gendir/\$( ls $gendir | grep \$genome | grep "\\.gff\$\\|\\.gff3\$" )
		check_gff.sh \$gff >> first_gff_checking.txt
	done
	# If the file gff_to_correct.txt is not empty, correct the GFF files.
	if [ -s gff_to_correct.txt ]
	then
		echo "WARNING : some GFF files are not formated as expected."
		echo "Trying to correct them..."
		mv $gendir ${gendir}_INIT
		cp -rL ${gendir}_INIT $gendir
		cat gff_to_correct.txt | while read gff
		do
			better_mRNA.py -i \$gff -o \$gff
		done
		echo "Done."

		cat valid_genomes.txt | while read genome
		do
			gff=$gendir/\$( ls $gendir | grep \$genome | grep "\\.gff\$\\|\\.gff3\$" )
			check_gff.sh \$gff >> gff_checking.txt
		done
	fi
	# If one or several GFF files are not correct :
	if [ -s gff_checking.txt ]
		then
		if grep -q ERROR gff_checking.txt
		then
			echo "ERROR : some GFF files are not formated as expected (see 'gff_checking.txt') :"
			cat gff_checking.txt
			exit 1
		fi
	fi


	# 3. Check the provided tree :
	if [ -s $tree ]
	then
		# Check that every genome in the list is in the phylogenetic tree. If not exit.
		# Also do the reverse and exit with respect to the prune_extra_taxa argument.
		# If set to 'True', a corrected tree will be generated, with only taxa also contained in the list. Otherwise, exit. 
		# (it can be convenient to download a tree with more genome than those of immediat interest).
		echo ""
		check_list_n_tree.py -names valid_genomes.txt -tree $tree -prune_extra_taxa True -out valid_tree.nwk
	else
		echo "ERROR : $tree DOES NOT EXISTS or IS EMPTY."
		exit 1
	fi


	# 4. Get phylip compatible names.
	# Get a two columns csv file with original species names as first col and PHYLIP compatible names as the second col.
	awk '{short=substr(\$0,1,5)NR; lim_short=substr(short,1,10); final=gensub(/_/,"","g",lim_short); print \$0"\t"final }' valid_genomes.txt > phylip_names.tsv
	# Get a new tree with PHYLIP compatible names.
	echo "Phylip tree :"
	get_phylip_named_tree.py -names phylip_names.tsv -tree valid_tree.nwk -out phylip_names_tree.nwk
	# In the TRG_table, use the phylip names instead of the original names.
	awk 'BEGIN{FS=OFS="\t"} FNR==NR {phylip_name[\$1]=\$2} FNR!=NR {if(FNR>1){\$1=phylip_name[\$1] ; \$3=phylip_name[\$3]} ; print \$0}' phylip_names.tsv $TRG_table > ${TRG_table}.tmp ; mv ${TRG_table}.tmp $TRG_table



	# 5. Built a channel for the rest of the pipeline
	cat phylip_names.tsv | while read line
	do
		genome=\$(echo "\${line}" | awk -F"\t" '{print \$1}')
		phy_genome=\$(echo "\${line}" | awk -F"\t" '{print \$2}')

		fasta=$gendir/\$( ls $gendir | grep \$genome | grep "\\.fas\$\\|\\.fna\$\\|\\.fasta\$\\|\\.fa\$" )
		gff=$gendir/\$( ls $gendir | grep \$genome | grep "\\.gff\$\\|\\.gff3\$" )

		# Creat a symlink to the fasta and gff files in the current directory.
		ln -s \$fasta \$phy_genome.fna
		ln -s \$gff \$phy_genome.gff

		echo "\${PWD}/\$phy_genome.fna__,__\${PWD}/\$phy_genome.gff" >> genome_files.txt
	done
	"""
}


process EXTRACT_CDS {

	input:
		tuple val(name), path(fasta), path(gff)
	output:
		tuple val(name), path(fasta), path(gff), path("${fasta}.fai"), path("${name}_CDS.fna"), path("${name}_CDS.faa")

	"""
	chmod -R +x ${projectDir}/bin

	# Use GffRead.

	# Need a fai index
	samtools faidx $fasta

	# Keep GFF sequences that are also in the FASTA file.
	awk 'BEGIN{FS=OFS="\t"} {print "^"\$1,""} END {print "^#"}' ${fasta}.fai | grep -f - $gff > gff_filterA

	# Also remove mRNA with undefined strand and features with abnormal end
	awk -F"\t" '
			FNR==NR {max[\$1]=\$2}
			FNR!= NR && ( /^#/ || (\$7 !~ /?/ && \$5 <= max[\$1]) )
			' ${fasta}.fai gff_filterA > gff_filterB

	# Extract the genomic CDS.
	# -V discard any mRNAs with CDS having in-frame stop codons
	# -x : write a fasta file with spliced CDS for each GFF transcript
	gffread -V -g $fasta -x ${name}_CDS.fna gff_filterB
	
	# Remove CDS missing a terminal stop codon and get a translated version of the FASTA.
	discard_CDS_missing_terminal_stop_codon.sh ${name}_CDS.fna
	"""
}


process ELONGATE_CDS {

	input:
		tuple val(name), path(gfasta), path(gff), path(fai), path(CDS_fna, stageAs: "CDS.fna"), path(CDS_faa, stageAs: "CDS.faa")
		
	output:
		path "${name}_CDS_elongated.fna", emit : CDS_elongated_fna
		path gfasta, emit : gfastas
		path fai, emit : fais
		
	"""
	chmod -R +x ${projectDir}/bin

	## In order to search for homologous CDS in the genome, get a FASTA with the elongated CDS of the genome (100 nucl upstream and downstream).
	# In case of multi matchs, the one with the best evalue is more likely to be the right one with this elongated version (bring genomic context).
	
	# Generate the elongated CDS FASTA.
	extend_cds.py \
	--gfasta $gfasta \
	--fai $fai \
	--gff $gff \
	--cds $CDS_fna \
	--size 99 \
	--mode nucl \
	--out ${name}_CDS_elongated.fna

	# Generate the elongated translated CDS FASTA.
	extend_cds.py \
	--gfasta $gfasta \
	--fai $fai \
	--gff $gff \
	--cds $CDS_fna \
	--size 99 \
	--mode simple \
	--out ${name}_CDS_elongated.faa
	"""
}


process ALIGNMENT_FASTA {

	input:
		path query_list
		path TRG_table
		val  focal
		path CDS_elongated
		path gfastas
		path fais

	output:
		// These output are optional since CDS with less than 2 neighbors with an homolog are filtered out at this step.
		path "*_toalign.fna", emit : fna, optional: true
		path "*.nwk", emit : newick, optional: true
		path "*_names_n_types.txt", emit : names_n_types, optional: true
		
	"""
	ls | grep "_CDS_elongated.fna\$" | sed "s/_CDS_elongated.fna\$//" > names.txt

	get_alignment_fasta.py $query_list $TRG_table names.txt $focal
	"""
}


process ANCESTRAL_ORFS {

	errorStrategy 'retry'

	publishDir "${params.outdir}/ssearch36_alignments", pattern: "*_ssearch36.tsv"
    publishDir "${params.outdir}/raw_ancorfs_fasta", pattern: "*_anc_ORFs.f*"

	input:
		path toali
		path tree
		path focal_CDS_faa

	output:
        path "*_anc_ORFs.f*", emit : raw_ancorfs_fastas
		path "*_ssearch36.tsv", emit : ssearch36_tsvs
		
	"""
	seq=\$(echo $toali | sed "s/_toalign.fna//")

	# Reconstruct the ancestral sequence.
	# <(awk '/^>/ {\$0=\$1} ; //' $toali) is just $toali with short headers.
	prank -d=<(awk '/^>/ {\$0=\$1} ; //' $toali) -t=$tree -o=\${seq} -once -showanc -prunetree

	# Select the node that corresponds to the most recent common ancestor of the focal genome and its closest outgroup neighbor(s).
	select_ancestor_node.py $tree $toali \${seq}.best.anc.dnd \${seq}.best.anc.fas -o \${seq}_anc.fna.tmp
	# Remove line breaks and '-' characters from the ancestral sequence.
	awk 'BEGIN {n=0;} /^>/ {if(n>0) printf("\\n%s\\n",\$0); else printf("%s\\n",\$0); n++; next;} {printf("%s",\$0);} END {printf("\\n");}' \${seq}_anc.fna.tmp | tr -d '-' > \${seq}_anc.fna ; rm \${seq}_anc.fna.tmp

	# Retrieve all the ORFs from the ancestral sequence.
	# --strand s (=sense), --type both (nucl and AA), --min_size 60 (minimal length 60 nucl)
	fasta_get_ORFs.py \${seq}_anc.fna --min_size 60 --strand s --type both --output \${seq}_anc_ORFs

	# Perform a pairwise alignment of each of the ancestral ORFs with the query CDS
	original_seq=\$(echo \${seq} | sed -e "s~__COLON__~:~g" -e "s~__SLASH__~/~g" -e "s~__EXCLAMATION__~!~g" -e "s~__PIPE__~|~g")
	faOneRecord $focal_CDS_faa \${original_seq} > \${seq}.faa
	# Ssearch with a tabular blast output format and and minimal e-value of 0.01
	ssearch36 -m8C -E 0.01 \${seq}.faa \${seq}_anc_ORFs.faa > \${seq}_anc_ORFs_vs_\${seq}_ssearch36.tsv
	"""
}


process ANCESTRAL_ORFS_MACSE {

	publishDir "${params.outdir}/ssearch36_alignments", pattern: "*_ssearch36.tsv"
    publishDir "${params.outdir}/raw_ancorfs_fasta", pattern: "*_anc_ORFs.f*"

	input:
		path toali
		path tree
		path focal_CDS_faa

	output:
        path "*_anc_ORFs.f*", emit : raw_ancorfs_fastas
		path "*_ssearch36.tsv", emit : ssearch36_tsvs
		
	"""
	seq=\$(echo $toali | sed "s/_toalign.fna//")

	# Align
	java -jar ${projectDir}/bin/macse_v2.07.jar -prog alignSequences -seq $toali -out_NT \${seq}_aligned.fna
	
	# Correct the NT file by replacing '!' characters by '-'
	sed 's|!|-|g' \${seq}_aligned.fna > \${seq}_aligned.fna.tmp
	# Remove description from the headers
	awk '/^>/ {print ">"substr(\$1, 2)} !/^>/ {print \$0}' \${seq}_aligned.fna.tmp > \${seq}_aligned.fna ; rm \${seq}_aligned.fna.tmp

	# Reconstruct the ancestral sequence.
	prank -d=\${seq}_aligned.fna -t=$tree -o=\${seq} -once -showanc -prunetree -keep

	# Select the node that corresponds to the most recent common ancestor of the focal genome and its closest outgroup neighbor(s).
	select_ancestor_node.py $tree $toali \${seq}.anc.dnd \${seq}.anc.fas -o \${seq}_anc.fna.tmp
	# Remove line breaks and '-' characters from the ancestral sequence.
	awk 'BEGIN {n=0;} /^>/ {if(n>0) printf("\\n%s\\n",\$0); else printf("%s\\n",\$0); n++; next;} {printf("%s",\$0);} END {printf("\\n");}' \${seq}_anc.fna.tmp | tr -d '-' > \${seq}_anc.fna ; rm \${seq}_anc.fna.tmp

	# Retrieve all the ORFs from the ancestral sequence.
	# --strand s (=sense), --type both (nucl and AA), --min_size 60 (minimal length 60 nucl)
	fasta_get_ORFs.py \${seq}_anc.fna --min_size 60 --strand s --type both --output \${seq}_anc_ORFs

	# Perform a pairwise alignment of each of the ancestral ORFs with the query CDS
	original_seq=\$(echo \${seq} | sed -e "s~__COLON__~:~g" -e "s~__SLASH__~/~g" -e "s~__EXCLAMATION__~!~g" -e "s~__PIPE__~|~g")
	faOneRecord $focal_CDS_faa \${original_seq} > \${seq}.faa
	# Ssearch with a tabular blast output format and and minimal e-value of 0.01
	ssearch36 -m8C -E 0.01 \${seq}.faa \${seq}_anc_ORFs.faa > \${seq}_anc_ORFs_vs_\${seq}_ssearch36.tsv
	"""
}


process ANCESTRAL_ORFS_MACSE_PHYML {

	// errorStrategy { task.attempt <= 1 ? 'retry' : 'ignore' }

	publishDir "${params.outdir}/ssearch36_alignments", pattern: "*_ssearch36.tsv"
    publishDir "${params.outdir}/raw_ancorfs_fasta", pattern: "*_anc_ORFs.f*"

	input:
		path toali
		path tree
		path focal_CDS_faa

	output:
        path "*_anc_ORFs.f*", emit : raw_ancorfs_fastas, optional: true
		path "*_ssearch36.tsv", emit : ssearch36_tsvs, optional: true
		
	"""
	seq=\$(echo $toali | sed "s/_toalign.fna//")

	# Align
	java -jar ${projectDir}/bin/macse_v2.07.jar -prog alignSequences -seq $toali -out_NT \${seq}_aligned.fna
	
	# Correct the NT file by replacing '!' characters by '-'
	sed 's|!|-|g' \${seq}_aligned.fna > \${seq}_aligned.fna.tmp
	# Remove description from the headers
	awk '/^>/ {print ">"substr(\$1, 2)} !/^>/ {print \$0}' \${seq}_aligned.fna.tmp > \${seq}_aligned.fna ; rm \${seq}_aligned.fna.tmp

	# Convert the NT alignment to a phylip file.
	fasta_to_phylip.py \${seq}_aligned.fna \${seq}_aligned.phylip

	# Get a subtree with only the genomes of the alignment.
	awk '/^>/ {print substr(\$1, 2)}' $toali > genome_names.txt
	tree_prune.py $tree genome_names.txt > subtree.nwk

	# Build a tree based on the alignment.
	### WARNING phyml needs a relative path after -i ###
	### WARNING depending on the computer used to install phyml 
	### (e.g. using nextflow with a singularity image built with another machine), 
	### launching phyml on another machine (or node, for HPCs) might prompt an "Illegal instruction" error, 
	### because CPUs are not configured the same way. 
	### See https://wiki.parabola.nu/Fixing_illegal_instruction_issues ###
	phyml -i \${seq}_aligned.phylip -d nt -v e -o lr -c 4 -a e -b 0 -f e -u subtree.nwk
	mv \${seq}_aligned.phylip_phyml_tree.txt alignment_tree.nwk

	# Reconstruct the ancestral sequence.
	prank -d=\${seq}_aligned.fna -t=alignment_tree.nwk -o=\${seq} -once -showanc -keep

	# Select the node that corresponds to the most recent common ancestor of the focal genome and its closest outgroup neighbor(s).
	select_ancestor_node.py $tree $toali \${seq}.anc.dnd \${seq}.anc.fas -o \${seq}_anc.fna.tmp
	# Remove line breaks and '-' characters from the ancestral sequence.
	awk 'BEGIN {n=0;} /^>/ {if(n>0) printf("\\n%s\\n",\$0); else printf("%s\\n",\$0); n++; next;} {printf("%s",\$0);} END {printf("\\n");}' \${seq}_anc.fna.tmp | tr -d '-' > \${seq}_anc.fna ; rm \${seq}_anc.fna.tmp

	# Check that the alignment tree do not imply to reconstruct an abberrant ancestor (i.e CDS matches in both bifurcations, the one with the focal, and the other).
	focal=\$(echo $focal_CDS_faa | sed "s/_CDS.faa//")
	ancestor_node=\$(head -n 1 \${seq}_anc.fna | sed "s/^>//")
	CDS_leaves=\$(grep "CDS" $toali | awk '{print substr(\$1, 2)}' | tr '\n' '|')
	CDS_leaves=\${CDS_leaves:-.}
	echo -e "orf\tnewick_string\tanode\tCDS_leaves" > input_check_alignment_tree.tsv
	echo -e "\${seq}\t\$(cat alignment_tree.nwk)\t\${ancestor_node}\t\${CDS_leaves}" >> input_check_alignment_tree.tsv
	check_alignment_tree.py input_check_alignment_tree.tsv output_check_alignment_tree.txt --focal_leaf_name \${focal}
	
	result_check_alignment_tree=\$(awk -F"\t" -v seq="\${seq}" '\$1 == seq {print \$2}' output_check_alignment_tree.txt)

	if [ "\$result_check_alignment_tree" == "FALSE" ]
	then
		echo "Stopping for \${seq} : the alignment tree implies to reconstruct an abberrant ancestor."
		exit 0
	fi

	# Retrieve all the ORFs from the ancestral sequence.
	# --strand s (=sense), --type both (nucl and AA), --min_size 60 (minimal length 60 nucl)
	fasta_get_ORFs.py \${seq}_anc.fna --min_size 60 --strand s --type both --output \${seq}_anc_ORFs

	# Perform a pairwise alignment of each of the ancestral ORFs with the query CDS
	original_seq=\$(echo \${seq} | sed -e "s~__COLON__~:~g" -e "s~__SLASH__~/~g" -e "s~__EXCLAMATION__~!~g" -e "s~__PIPE__~|~g")
	faOneRecord $focal_CDS_faa \${original_seq} > \${seq}.faa
	# Ssearch with a tabular blast output format and and minimal e-value of 0.01
	ssearch36 -m8C -E 0.01 \${seq}.faa \${seq}_anc_ORFs.faa > \${seq}_anc_ORFs_vs_\${seq}_ssearch36.tsv
	"""
}


process ANCESTRAL_ORFS_RAXML {

	publishDir "${params.outdir}/ssearch36_alignments", pattern: "*_ssearch36.tsv"
    publishDir "${params.outdir}/raw_ancorfs_fasta", pattern: "*_anc_ORFs.f*"

	input:
		path toali
		path tree
		path focal_CDS_faa

	output:
        path "*_anc_ORFs.f*", emit : raw_ancorfs_fastas
		path "*_ssearch36.tsv", emit : ssearch36_tsvs
		
	"""
	seq=\$(echo $toali | sed "s/_toalign.fna//")

	# Align
	java -jar ${projectDir}/bin/macse_v2.07.jar -prog alignSequences -seq $toali -out_NT \${seq}_aligned.fna
	
	# Correct the NT file by replacing '!' characters by '-'
	sed 's|!|-|g' \${seq}_aligned.fna > \${seq}_aligned.fna.tmp
	# Remove description from the headers
	awk '/^>/ {print ">"substr(\$1, 2)} !/^>/ {print \$0}' \${seq}_aligned.fna.tmp > \${seq}_aligned.fna ; rm \${seq}_aligned.fna.tmp

	# Build a tree based on the alignment, while preserving the original whole-genome tree topology.
	tree_prune_with_fasta.py $tree \${seq}_aligned.fna ${tree}_pruned
	raxmlHPC -s \${seq}_aligned.fna -m GTRGAMMA -g ${tree}_pruned -n \${seq} -p 12345

	# Reconstruct the ancestral sequence.
	prank -d=\${seq}_aligned.fna -t=RAxML_bestTree.\${seq} -o=\${seq} -once -showanc -keep

	# Select the node that corresponds to the most recent common ancestor of the focal genome and its closest outgroup neighbor(s).
	select_ancestor_node.py $tree $toali \${seq}.anc.dnd \${seq}.anc.fas -o \${seq}_anc.fna.tmp
	# Remove line breaks and '-' characters from the ancestral sequence.
	awk 'BEGIN {n=0;} /^>/ {if(n>0) printf("\\n%s\\n",\$0); else printf("%s\\n",\$0); n++; next;} {printf("%s",\$0);} END {printf("\\n");}' \${seq}_anc.fna.tmp | tr -d '-' > \${seq}_anc.fna ; rm \${seq}_anc.fna.tmp

	# Retrieve all the ORFs from the ancestral sequence.
	# --strand s (=sense), --type both (nucl and AA), --min_size 60 (minimal length 60 nucl)
	fasta_get_ORFs.py \${seq}_anc.fna --min_size 60 --strand s --type both --output \${seq}_anc_ORFs

	# Perform a pairwise alignment of each of the ancestral ORFs with the query CDS
	original_seq=\$(echo \${seq} | sed -e "s~__COLON__~:~g" -e "s~__SLASH__~/~g" -e "s~__EXCLAMATION__~!~g" -e "s~__PIPE__~|~g")
	faOneRecord $focal_CDS_faa \${original_seq} > \${seq}.faa
	# Ssearch with a tabular blast output format and and minimal e-value of 0.01
	ssearch36 -m8C -E 0.01 \${seq}.faa \${seq}_anc_ORFs.faa > \${seq}_anc_ORFs_vs_\${seq}_ssearch36.tsv
	"""
}


process ANCESTRAL_ORFS_PREQUEL {

	publishDir "${params.outdir}/ssearch36_alignments", pattern: "*_ssearch36.tsv"
    publishDir "${params.outdir}/raw_ancorfs_fasta", pattern: "*_anc_ORFs.f*"

	input:
		path toali
		path tree
		path focal_CDS_faa

	output:
        path "*_anc_ORFs.f*", emit : raw_ancorfs_fastas
		path "*_ssearch36.tsv", emit : ssearch36_tsvs
		
	"""
	seq=\$(echo $toali | sed "s/_toalign.fna//")

	# Align
	java -jar ${projectDir}/bin/macse_v2.07.jar -prog alignSequences -seq $toali -out_NT \${seq}_aligned.fna
	
	# Correct the NT file by replacing '!' characters by '-'
	sed 's|!|-|g' \${seq}_aligned.fna > \${seq}_aligned.fna.tmp
	# Remove description from the headers
	awk '/^>/ {print ">"substr(\$1, 2)} !/^>/ {print \$0}' \${seq}_aligned.fna.tmp > \${seq}_aligned.fna ; rm \${seq}_aligned.fna.tmp

	# Reconstruct the ancestral sequence.
	phyloFit --seed 12546582 --msa \${seq}_aligned.fna --tree phylip_names_tree.nwk
	tree_prune_with_fasta.py phylip_names_tree.nwk \${seq}_aligned.fna pruned_tree.nwk
	gotree_amd64_linux rename -i pruned_tree.nwk -l 1 --tips=false --internal -a | sed 's/N000/N/g' > final_tree.nwk
	cat  phyloFit.mod | sed "s/^TREE\:.*/TREE: $(cat final_tree.nwk)/" > phylofit_corTree.mod
	prequel \${seq}_aligned.fna phylofit_corTree.mod \${seq}
	prequel \${seq}_aligned.fna phylofit_corTree.mod \${seq} -n
	
	# Select the node that corresponds to the most recent common ancestor of the focal genome and its closest outgroup neighbor(s).
	select_ancestor_node.py $tree $toali final_tree.nwk \${seq}.anc.fas -o \${seq}_anc.fna.tmp
	# Remove line breaks and '-' characters from the ancestral sequence.
	awk 'BEGIN {n=0;} /^>/ {if(n>0) printf("\\n%s\\n",\$0); else printf("%s\\n",\$0); n++; next;} {printf("%s",\$0);} END {printf("\\n");}' \${seq}_anc.fna.tmp | tr -d '-' > \${seq}_anc.fna ; rm \${seq}_anc.fna.tmp

	# Retrieve all the ORFs from the ancestral sequence.
	# --strand s (=sense), --type both (nucl and AA), --min_size 60 (minimal length 60 nucl)
	fasta_get_ORFs.py \${seq}_anc.fna --min_size 60 --strand s --type both --output \${seq}_anc_ORFs

	# Perform a pairwise alignment of each of the ancestral ORFs with the query CDS
	original_seq=\$(echo \${seq} | sed -e "s~__COLON__~:~g" -e "s~__SLASH__~/~g" -e "s~__EXCLAMATION__~!~g" -e "s~__PIPE__~|~g")
	faOneRecord $focal_CDS_faa \${original_seq} > \${seq}.faa
	# Ssearch with a tabular blast output format and and minimal e-value of 0.01
	ssearch36 -m8C -E 0.01 \${seq}.faa \${seq}_anc_ORFs.faa > \${seq}_anc_ORFs_vs_\${seq}_ssearch36.tsv
	"""
}


process ANCORFS_FASTA {

    publishDir "${params.outdir}/ancorfs_fasta/", mode : 'copy'

	input:
		path raw_ancorfs_fastas
		path ssearch36_tsvs

	output:
		path "*"
		
	"""
	get_ancorfs_fasta.py $ssearch36_tsvs --mode best --mode 1e-3 --mode 1e-2
	"""
}
