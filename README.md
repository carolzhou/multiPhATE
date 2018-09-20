/Code3/ - multiPhATE (beta version)

# This code was developed by Carol L. Ecale Zhou and Jeffrey Kimbrel at Lawrence Livermore National Laboratory.
# THIS CODE IS COVERED BY THE BSD LICENSE. SEE INCLUDED FILE BSD-3.pdf FOR DETAILS.

ABOUT THE MULTI-PHATE PIPELINE

multiPhATE is a throughput version of PhATE, which is described below. The multiPhate.py code takes a single argument (myMultiPhate.config; use multiPhate.config as an example) and uses it to generate a phate.config file for each genome being annotated. Then, multiPhate.py invokes the PhATE pipeline (via phate_runPipeline.py) for each genome listed in the multiPhate.config file.

ABOUT THE PHATE PIPELINE

PhATE is a fully automated computational pipeline for annotating phage genomes or for identifying putative phage genes in bacterial genomes. PhATE is written in Python 2.7. Code execution is controled by a configuration file (see sample.config), which can be tailored to run specific gene finders and to blast sequences against specific phage- virus-centric data sets, in addition to more generic (genome, protein) data sets. PhATE runs at least one gene finding algorithm, then annotates the genome, gene, and protein sequences using blast and a set of fasta sequence databases, and uses an hmm search against the pVOG database. If more than one gene finder is run, PhATE will also provide a side-by-side comparison of the genes called by each gene caller. The user specifies the preferred gene caller, and the genes and proteins predicted by that caller are annotated using blast against the supporting databases. Classification of each protein sequence into a pVOG group is followed by generation of an alignment-ready fasta file.

HOW TO SET UP PHATE ON YOUR LOCAL MACHINE

Place multiPhate.py and phate_runPipeline.py in your main execution directory. Create two subdirectories:  PipelineInput/ and PipelineOutput/
Place your phage genome fasta files (myGenome.fasta) into the PipelineInput/ subdirectory (single contig)
Place your configuration file (multiPhate.config) in the main directory (same level as multiPhate.py)

You will need to acquire one or more of the databases listed below under SUPPORING DATABASES, and the 3rd party codes listed under SUPPORTING 3rd PARTY CODES. You will need to acquire at least one of the supported gene finders, but it is recommended to run as many of the four gene finders as is feasible so that the results can be more meaningfully compared. You will need to specifiy the locations of the supporting data sets and codes in phate_runPipeline.py under "USER CONFIGURATION", and you will need to locate your genome file to the PipelineInput directory. Place your configuration file at the top level (alongside phate_runPipeline.py).

You may up- or down-regulate verbosity in the phate_runPipeline.py file, under "USER CONFIGURATION". This includes an option to clean the (voluminous) raw blast and hmm search data from the output directories.

(You will need to modify the sample configuration file, multiPhate.config, accordingly.)
Run the PhATE pipeline at the command line using the following command:
$ python multiPhate.py myMultiPhate.config  

SUPPORTING DATABASES

It is recommended that the user acquire as many of the following sequence databases and associated codes as is feasible, although none are actually required to run the code. (You may specify "translate_only='true'" to do gene finding and stop at that point.) Databases are listed with at least one way of acquiring them, but there may be additional sources, and it is possible to substitute subsets of blast databases (e.g., a subset of the NCBI gene database inplace of Refseq Gene).  

* NCBI virus genomes - ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/ or https://www.ncbi.nlm.nih.gov/genome/viruses/viral.1.1.1.genomic.fna.gz
* NCBI accession2taxid file - ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/nucl_gb.accession2taxid.gz or https://www.ncbi.nlm.nih.gov/guide/taxonomy/nucl_gb.accession2taxid.gz
* NCBI Refseq Protein - download using blast+: /bin/update_blastdb.pl refseq_protein
* NCBI Refseq Gene - download using blast+: /bin/update_blastdb.pl refseqgene. This database contains primarily human sequences. To acquire bacterial gene sequences you may need to download them via http://ncbi.nlm.nih.gov/gene and process the data to generate a fasta data set.
* NCBI Swissprot - download using blast+: /bin/update_blastdb.pl swissprot 
* NR - ftp://ftp.ncbi.nlm.nih.gov/nr/
* KEGG virus subset - (available by license) http://www.kegg.jp/kegg/download/
* KEGG associated files - T40000.pep, T40000.nuc, vg_enzyme.list, vg_genome.list, vg_ko.list, vg_ncbi-geneid.list, vg_ncbi-proteinid.list, vg_pfam.list, vg_rs.list, vg_tax.list, vg_uniprot.list 
* Phantome protein fasta sequences - http://www.phantome.org/Downloads/phage_proteins_nnnnnnnnn.fasta
* pVOGs prepared database (pVOGs.faa) - included in PhATE distribution. This data set was derived by C. Zhou from the pVOGs fasta database. For use in PhATE, the sequence fasta headers have been modified to include the pVOG identifiers (all groups to which each sequence belongs). This re-formatting facilitates pVOG group identification and construction of the alignment-ready fasta files. Codes for reconstructing this modified data set are included in the PhATE distribution. Note that the pVOGs are not mutually exclusive, meaning that a sequence may have membership in more than one VOG group. The codes included in the phate distribution will combine identifiers that belong to a given sequence and list all the VOG identifiers in the fasta header. In this way, the pVOG fasta database used in PhATE is also non-redundant. See documentation in DatabasePrep/dbPrep_createPvogFastaFile.py for instructions how to update your local pVOGs data set for use in PhATE, but you can start with the pVOGs.faa file included in the PhATE distribution. 

For simplicity in configuring the locations of dependent databases in the phate_runPipeline.py driver script, it is suggested that the above databases be placed in a directory structure as follows:
/Databases/
	KEGG/
	NCBI/
		Virus_Genome/
		Virus_Protein/
	NR/
	Phantome/
	Refseq/
		Protein/
		Gene/
	Swissprot/
	pVOGs/

SUPPORTING 3rd PARTY CODES

* BioPython - https://biopython.org/wiki/Download 
* EMBOSS package - https://sourceforge.net/directory/os:mac/?q=EMBOSS
* Blast+ https://ncbi.nlm.nih.gov/blast
* GeneMarkS - http://exon.gatech.edu/Genemark/index.html (optional; available by license)
* Glimmer - https://ccb.jhu.edu/software/glimmer/ (optional)
* Prodigal - https://github.com/hyattpd/Prodigal (optional)
* PHANOTATE - https://github.com/deprekate/PHANOTATE
* jackhmmer - https://www.ebi.ac.uk/Tools/hmmer/search/jackhmmer

HOW TO WRITE A CONFIGURATION FILE

Locations of supporting databases and codes are to be specified in the file, multiPhate.py, according to their locations on your local system. 

A sample configuration file is provided (multiPhate.config). Make a copy of this file and rename it accordingly (eg., myMultiPhate.config). The sample.config file is configured according to established default parameters (just about everything turned 'off' initially). Any of the parameters may be modified (switches turned 'on' or 'off') by assigning 'true' or 'false'. It is suggested that you turn swithes off, then install each supporting database in turn and test the code. A switch not included in the configuration file is CLEAN_RAW_DATA, which you will find in multiPhate.py (under USER CONFIGURATION). By default all raw blast and hmm outputs (written under PipelineOutput/BLAST/ and PipelineOutput/HMM/ subdirectories) will be deleted, so as to limit output volume. However, if you wish to retain that output, turn the CLEAN_RAW_DATA switch to 'True' in phate_runPipeline.py.

PHATE PIPELINE OUTPUT FILES
In the user-specified output directory, under PipelineOutput/, the following files will be written:
* gene.fnt and protein.faa fasta files generated using the designated preferred gene finder
* CGC_results.txt - side-by-side comparison of all gene finder results (if at least two were run)
* phate_sequenceAnnotation_main.out - tabbed integrated annotation results (also written in gff format)
* BLAST/ directory - raw blast results for genome (under Genome/) and proteins (under Protein/) (if CLEAN_RAW_DATA switch in phate_runPipeline.py is 'False')
* pVOG groupings comprising alignment-ready fasta files (under BLAST/Protein/) based on best blast hits
* HMM/ directory - raw hmm search results for protein sequences (under Protein/) (if CLEAN_RAW_DATA switch in phate_runPipeline.py is 'False')
* pVOG grouping comprising alignment-ready fasta files (under HMM/Protein/) based on best hmm hits
* log files capturing details of the processing and time stamps
* your myPhate.config file, to record exactly how you configured the pipeline for the current run.

RUNNING PHATE AS AN "EMBARASSINGLY PARALLEL" CODE

Pipeline outputs are written to user-specified output directories (specified in your myMultiPhate.config file, and in the auto-generated myPhate.config files, one for each genome). In this way, PhATE may be set up to run in parallel on any number of nodes of a compute cluster without fear of clashes in writing results. Implementing a code in parallel is system dependent.

FURTHER RECOMMENDATIONS

PhATE was originally developed for anotating phage genome sequences. However, PhATE may also be useful in helping to identify phage genes within bacterial genomes (i.e., prophage). Thus, the user has the option of running multiple gene callers for bacterial genome sequence (GeneMarkS, Glimmer, Prodigal) and a new gene finder specifically for phage (PHANOTATE). The calls from any two or more of these callers are compared by the PhATE/CGC code so that the user can examine the calls that agree or disagree among the callers, and then run PhATE again, selecting the caller of choice for annotating the sequence. 

Although most supporting data sets are phage- or virus-centric, NR and Refseq Protein are included in order to help identify genes/functions that are not known in the virus/phage gene data sets. However, PhATE is not intended to be the sole source for anntotation of bacterial genome sequences, as PhATE is tailored for identification of genes and functions in phage.

The NR database has grown enormously large. It is recommended to use a smaller database, such as Refseq Protein in place of NR. Furthermore, annotating with NR will add greatly to the time required for processing a genome through PhATE. It is recommended that NR be turned 'off' ('false') until one desires to preform a full/final annotation of the genome of interest, if using NR.

Because the behavior of 3rd party codes can sometimes be unpredictable, it is recommended that the user replace spaces and unusual characters in their fasta headers with the underscore character.

PLANNED FURTHER DEVELOPMENT

Plans include adding additional functionality to the pipeline, such as HMMER and additional hmm databases.
# multiPhATE
