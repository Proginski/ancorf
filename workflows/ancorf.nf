/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WELCOME
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

log.info "\n"
log.info "Welcome to DENSE."
log.info "\n"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PARAMETERS MANAGMENT
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { validateParameters; paramsHelp; paramsSummaryLog; fromSamplesheet } from 'plugin/nf-validation'

// Print help message, supply typical command line usage for the pipeline
if (params.help) {
   log.info paramsHelp("nextflow run proginski/ancorf --gendir <DIR WITH GFF AND FASTA> --focal <FOCAL_GENOME_NAME> --tree <NEWICK WITH FOCAL AND NEIGHBORS> --outdir <OUTDIR>")
   exit 0
}

// Validate input parameters
validateParameters()

// Print summary of supplied parameters
log.info paramsSummaryLog(workflow)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CHECK_INPUTS               } from '../modules/local/ancorf_modules.nf'
include { EXTRACT_CDS                } from '../modules/local/ancorf_modules.nf'
include { ELONGATE_CDS               } from '../modules/local/ancorf_modules.nf'
include { ALIGNMENT_FASTA            } from '../modules/local/ancorf_modules.nf'
include { ANCESTRAL_ORFS             } from '../modules/local/ancorf_modules.nf'
include { ANCESTRAL_ORFS_MACSE       } from '../modules/local/ancorf_modules.nf'
include { ANCESTRAL_ORFS_MACSE_PHYML } from '../modules/local/ancorf_modules.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


/*
This nextflow pipeline reconstructs the ancestral sequence of CDS and identifies the ORF(s) that gave birth to the present CDS.

In order to infer evolutionnary relationships, it requires for each CDS the presence of an homolog in at least two neighbor genomes (CDS or intergenic region).
As it is part of a wider analysis, the current pipeline assumes the existence of :
	- ${PWD}/GENOMES/${genome}.fna 
	- ${PWD}/GENOMES/${genome}.gff

Warnings : 
	- the input newick tree must not contain duplicated taxon labels
	- in ${PWD}/GENOMES/${genome}.gff : the CDS having as "Parent" a features different from mRNA might be discarded by gffread.
*/


// WORKFLOW
workflow ANCORF {

	/*
	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	Check the inputs before the beginning.
	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	*/

	CHECK_INPUTS(
				 file(params.gendir),
				 file(params.tree),
				 file(params.trg_table)
				)

	// genomic FASTA and GFF3 pairs that will be processed
	genome_ch = CHECK_INPUTS.out.genome_files
		.splitText() 
		.map { tuple( it.strip().split("__,__") ) }
		.map { fasta, gff -> [ file(fasta).getBaseName(), fasta, gff ] }

	/*
	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	Exctract the CDS of each genome
	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	*/

	EXTRACT_CDS(genome_ch)

	// Also get an elongated version of the trasnlated CDS for every genome (subjects).
	ELONGATE_CDS(EXTRACT_CDS.out)
	ELONGATE_CDS.out.CDS_elongated_fna.collect().set { elongate_cds_fastas }
	ELONGATE_CDS.out.gfastas.collect().set { elongate_cds_gfastas }
	ELONGATE_CDS.out.fais.collect().set { elongate_cds_fais }

	// Get the new, phylip compliant, names of the focal genome.
	focal = CHECK_INPUTS.out.phylip_names
		.splitText() 
		.map { tuple( it.strip().split("\t") ) }
		.filter { it[0] == params.focal }
		.map { it[1] }
		.first()

	// Get the name of the focal CDS aa FASTA file
	focal_CDS_faa = EXTRACT_CDS.out
		.filter { it[0] == focal.value }
		.map { it[5] }
		.first()

	/*
	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	Get a FASTA with the CDS and its neighbor nucleotide homologous sequences.
	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	*/

	ALIGNMENT_FASTA(
		file(params.queries), 
		CHECK_INPUTS.out.trg_table,
		focal,
		elongate_cds_fastas,
		elongate_cds_gfastas,
		elongate_cds_fais
		)

	/*
	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	Reconstruct the ancestral ORFs.
	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	*/
	
	if (params.mode == "genomes_not_aligned") {
		ANCESTRAL_ORFS( 
			ALIGNMENT_FASTA.out.fna.flatten(),
			CHECK_INPUTS.out.tree,
			focal_CDS_faa
			)
	} else if (params.mode == "genomes_aligned") {
		ANCESTRAL_ORFS_MACSE( 
			ALIGNMENT_FASTA.out.fna.flatten(),
			CHECK_INPUTS.out.tree,
			focal_CDS_faa
			)
	} else if (params.mode == "genes_aligned") {
		ANCESTRAL_ORFS_MACSE_PHYML( 
			ALIGNMENT_FASTA.out.fna.flatten(),
			CHECK_INPUTS.out.tree,
			focal_CDS_faa
			)
	}

 }
