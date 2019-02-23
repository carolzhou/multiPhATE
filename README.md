# multiPhATE v.0.5
/MultiPhate/ - multiPhATE (beta version)

This code was developed by Carol L. Ecale Zhou and Jeffrey Kimbrel at Lawrence Livermore National Laboratory.

THIS CODE IS COVERED BY THE BSD LICENSE. SEE INCLUDED FILE BSD-3.pdf FOR DETAILS.

ABOUT THE MULTI-PHATE PIPELINE DRIVER

MultiPhATE is a throughput version of PhATE, which is described below. The multiPhate.py code is a command-line program that takes a single argument (hereafter referred to as, "multiPhate.config"; use sample.multiPhate.config as starting point) and uses it to generate a phate.config file (suitably named) for each genome being annotated. Then, multiPhate.py invokes the PhATE pipeline (via phate_runPipeline.py) for each genome.

ABOUT THE PHATE PIPELINE

PhATE is a fully automated computational pipeline for identifying and annotating phage genes in genome sequence. PhATE is written in Python 3.7, and runs on Linux and Mac operating systems. Code execution is controled by a configuration file, which can be tailored to run specific gene finders and to blast sequences against specific phage- and virus-centric data sets, in addition to more generic (genome, protein) data sets. PhATE runs at least one gene finding algorithm, then annotates the genome, gene, and protein sequences using blast and a set of fasta sequence databases, and uses an hmm search against the pVOG database. If more than one gene finder is run, PhATE will also provide a side-by-side comparison of the genes called by each gene caller. The user specifies the preferred gene caller, and the genes and proteins predicted by that caller are annotated using blast against the supporting databases. Classification of each protein sequence into a pVOG group is followed by generation of an alignment-ready fasta file. By convention, genome sequence files end with extension, ".fasta"; gene nucleotide fasta files end with, ".fnt", and cds amino-acid fasta files end with, ".faa".

HOW TO SET UP MULTI-PHATE ON YOUR LOCAL MACHINE

First, create a working directory on your computer for running multiPhATE. Then, acquire the multiPhATE package from github. This can be done either by downloading a zip file directly from the multiPhATE repository, or by cloning the repository. The first method is recommended, but the second is certainly an option:

*) To download the zip file:  Use a browser and navigate to https://github.com/carolzhou/multiPhATE. Press the green button "Clone or download", and download the zip file. Then, unzip the package in your working (main execution "multiPhate") directory.

$ cd myMultiphateDir

$ unzip multiPhate-master.zip

*) To clone from github:  Acquire git from https://git-scm.com/downloads. Naviate to your working (main execution "multiPhATE") directory, and clone multiPhATE from the command line: 

$ git init

$ git clone https://github.com/carolzhou/mulitPhATE

(Complete instructions for using git and github can be found at http://help.github.com.)

Now, be sure that multiPhate.py and phate_runPipeline.py and associated files and directories are in your main execution "multiPhATE" directory. Check that the two subdirectories: PipelineInput/ and PipelineOutput/ are present (should already exist in the downloaded distribution). Place your phage genome fasta files (genome1.fasta, genome2.fasta, etc.) into the PipelineInput/ subdirectory. Place your configuration file (ie, your copy of sample.multiPhate.config) in the main execution directory (same level as multiPhate.py). A word of caution here:  it is always best to name your files and fasta contigs as strings lacking any spaces or special characters, as third-party codes over which we have no control may balk when encountering odd characters or spaces. 

You will need to acquire one or more of the databases listed below under SUPPORING DATABASES (Phantome and pVOGs are included in the multiPhATE distribution, so it is possible to begin with just those), and the 3rd party codes listed under SUPPORTING 3rd PARTY CODES. You will need to acquire at least one of the supported gene finders, but it is recommended to run as many of the four gene finders as is feasible so that the results can be more meaningfully compared. You will need to specifiy the locations of the supporting data sets and codes in the multiPhATE config file (see multiPhate.config), and you will need to locate your genome file(s) to the PipelineInput/ subdirectory. Once you have acquired the third-party codes and databases, you will be ready to configure the multiPhate.config file.

HOW TO WRITE A CONFIGURATION FILE

Summary:
Availability and locations of supporting databases and codes are to be specified in a configuration file. A sample configuration file is provided, called "sample.multiPhate.config". Make a copy of this file and rename it accordingly (eg., myGenomeSet_multiPhate.config). Hereafter we refer to this file as, multiPhate.config. The multiPhate.config file is configured according to established default parameters (just about everything turned off initially). Any of the parameters may be modified (switches turned on or off) by assigning 'true' or 'false'. It is suggested that you turn swithes off, then install each supporting gene finder and database in turn and test the pipeline.

Procedure:

1) At the command line, make a copy of the file, sample.multiPhate.config, and name it appropriately (hereafter referred to as 'multiPhate.config'):  $ cp sample.multiPhate.config multiPhate.config.  Then, edit your config file as described below.

2) List of Genomes:
For each genome to be processed, provide six lines under "Genome List:" and before "END of list":  for each genome, you need to list the genome number, the name of the genome fasta file, the genome type (typically 'phage', but could be 'bacteria'), the species, if known (no spaces), the name of the genome, and a name for the output directory to hold this genome's output files (again, no spaces), in that order. You can simply copy/paste the six lines provided as many times as needed, and fill in the information appropriate for each genome.

3) Processing Information:
You may configure the pipeline to perform gene finding only, or gene finding plus functional annotation. For example, you may want to examine the results of multiple gene finders before going forward with functional annotation. In order to configure phate to run gene finding only, set translate_only to 'true'; in this way, only gene-calling and translation (to peptide sequence) will be performed. If you set translate_only to 'false', then the pipeline will not stop at the translation step, but will proceed with functional annotation of the predicted genes (ie, blast and/or hmm). Normally the genetic_code should be set to '11', for prokaryotic.

4) Gene Callers:
The gene_caller option specifies which gene caller's results (ie, gene calls) will be used for subsequent functional annotation. The choices are:  'phanotate', 'genemarks', 'prodigal', or 'glimmer'.  To run a gene caller, you must have acquired that third-party code and installed it locally for use with multiPhATE. For each gene caller you wish to have run, set that caller's parameter to 'true'. In the usual case, you will want to specify gene_caller='phanotate' for annotation of phage genomes.

5) Annotation:
Set to 'true' each blast or hmm process that you want to be run. Note that you must have acquired the associated database, and in the next section (Databases) you must configure the location of each database. You may also set the desired blast parameters. The blast_identity sets the minimum identity that will be considered; any blast result below that threshold will be ignored. The hit_count parameters will determine how many top hits will be reported. Currently the only hmm_program that is supported by multiPhate is 'jackhmmer', and it is only run with the pVOGs database (future releases of multiPhate are expected to support additional hmm analyses).

6) Databases:
For each database that you have in-house, specify the full path/filename. Note that you must prepare in advance all blast databases by running the "makeblastdb" utility (see instructions with blast+ code for how to do that). MultiPhate will only run with blast+; it does not support legacy blast. For instructions where to download the databases, see the SUPPORTING DATABASES section below. Note that KEGG is available by license. Note also that in some cases additional files are required. In this case, place the additional file(s) in the same directory as the associated blast database. For example, place the NCBI accession2taxid file in the same directory as your NCBI virus genome file (see below). If you are downloading datasets that you anticipate using specifically with multiPhATE, then it is suggested, for convenience, that you save them in the Databases/ folder in the multiPhATE distribution, but any database can be located anywhere on your local system; you need only indicate in the multiPhate.config file the full path/filename for each database. Remember, the pVOGs and Phantome data sets are included in the multiPhATE distribution in the Databases/ folder, but you will need to run makeblastdb to render the datasets blast-able ($ makeblastdb -help).

7) Verbosity:
You may up- or down-regulate verbosity in the multiPhate.config file, under "# VERBOSITY". This includes an option to clean the (voluminous) raw blast and hmm search data from the output directories. It is suggested that clean_raw_data, phate_progress, and cgc_progress be set to 'true'. The warnings and messages, when set to 'true', will generate voluminous output; set these to 'true' when trouble-shooting the pipeline.


PIPELINE EXECUTION

Run the PhATE pipeline at the command line by passing your multiPhate.config file as an argument to the multiPhate.py pipeline driver script, as follows: $ python multiPhate.py multiPhate.config


SUPPORTING DATABASES

It is recommended that the user acquire as many of the following sequence databases and associated codes as is feasible, although none are actually required to run the code. (You may specify "translate_only='true'" to do gene finding then translation, and then stop at that point.) Databases are listed with at least one way of acquiring them, but there may be additional sources, and it is possible to substitute subsets of blast databases (e.g., a subset of the NCBI gene database in place of Refseq Gene).

In some cases, you will download a database through a web interface, and in other cases you may use blast+ to download a database at the command line, if the database of interest is provided by NCBI. The latter method may be more convenient for larger data sets (eg, NCBI Refseq Protein). Blast+ includes a script (/bin/update_blastdb.pl), which can be invoked at the command line to automatically download a specified database. In order for blast+ to search a database, the database must first be converted to a blast-able object using the blast+ program, makeblastdb. Once you have installed blast+, you can query for help in using the program. For example, type at the command line: $ makeblastdb -help.

NCBI virus genomes - ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/ or https://www.ncbi.nlm.nih.gov/genome/viruses/viral.1.1.1.genomic.fna.gz

NCBI-associated file:  accession2taxid file - ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/nucl_gb.accession2taxid.gz or https://www.ncbi.nlm.nih.gov/guide/taxonomy/nucl_gb.accession2taxid.gz

NCBI Refseq Protein - download using blast+: /bin/update_blastdb.pl refseq_protein

NCBI Refseq Gene - download using blast+: /bin/update_blastdb.pl refseqgene. This database contains primarily human sequences. To acquire bacterial gene sequences you may need to download them via http://ncbi.nlm.nih.gov/gene and process the data to generate a fasta data set. Support for doing this is not provided with the multiPhATE distribution.

NCBI Swissprot - download using blast+: /bin/update_blastdb.pl swissprot

NR - ftp://ftp.ncbi.nlm.nih.gov/nr/

KEGG virus subset - (available by license) http://www.kegg.jp/kegg/download/

KEGG associated files - T40000.pep, T40000.nuc, vg_enzyme.list, vg_genome.list, vg_ko.list, vg_ncbi-geneid.list, vg_ncbi-proteinid.list, vg_pfam.list, vg_rs.list, vg_tax.list, vg_uniprot.list

Phantome protein fasta sequences - http://www.phantome.org/Downloads/phage_proteins_nnnnnnnnn.fasta. (A version of Phantome is included in the multiPhATE distribution.)

pVOGs prepared database (pVOGs.faa) - included in PhATE distribution. This data set was derived by C. Zhou from the pVOGs fasta database. For use in PhATE, the sequence fasta headers have been modified to include the pVOG identifiers (all groups to which each sequence belongs). This re-formatting facilitates pVOG group identification and construction of the alignment-ready fasta files. Codes for reconstructing this modified data set are included in the PhATE distribution. Note that the pVOGs are not mutually exclusive, meaning that a sequence may have membership in more than one VOG group. The codes included in the phate distribution will combine identifiers that belong to a given sequence and list all the VOG identifiers in the fasta header. In this way, the pVOG fasta database used in PhATE is also non-redundant. See documentation in DatabasePrep/dbPrep_createPvogFastaFile.py for instructions how to update your local pVOGs data set for use in PhATE, but you can start with the pVOGs.faa file included in the PhATE distribution.

For simplicity in configuring the locations of dependent databases in the multiPhate.config file, it is suggested that the above databases be placed in a directory structure as follows: 

Databases/
 
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

You must specify in your multiPhate.config file the locations of the data sets that you will be using. Although it is recommended that you place your databases in the above directory structure, they can reside anywhere locally on disk, but in any case you must specify the full directory path/filename to a given resource in your multiPhate.config file.


SUPPORTING 3rd PARTY CODES

Note that some third-party codes are required for multiPhATE, but others are optional, as indicated below. Some of these codes can be installed in a Conda environment. Codes that can be installed via Conda are so indicated below. If using Conda, follow the instructions that occur at the bottom of this section. Otherwise, install these codes globally, following the instructions provided with each package from the source.

BioPython - https://biopython.org/wiki/Download (required) (conda)

EMBOSS package - https://sourceforge.net/directory/os:mac/?q=EMBOSS (required) (conda)

Blast+ https://ncbi.nlm.nih.gov/blast (optional) (conda)

GeneMarkS - http://exon.gatech.edu/Genemark/index.html (optional; available by license)

Glimmer - https://ccb.jhu.edu/software/glimmer/ Use Glimmer version 3. (optional) (conda)

Prodigal - https://github.com/hyattpd/Prodigal (optional) (conda)

PHANOTATE - A Python 3-compatible version of PHANOTATE is included in the multiPhATE distribution, in the ExternalCodes/ folder. Unzip PHANOTATE.zip, and follow instructions in the README. Future updates to PHANOTATE will be available at https://github.com/deprekate/PHANOTATE. (optional)

jackhmmer - https://www.eddylab.org/software.html or http://www.hmmer.org/download.html. Download HMMER; jackhmmer is included in this package (optional) (conda - hmmer) 

tRNAscan-SE - https://www.eddylab.org/software.html - select tRNAscan-SE download link (conda)

Third-party codes should be installed globally whenever possible. However, it is recommended that PHANOTATE be installed under the ExternalCodes/ subdirectory in the execution/working directory. (The ExternalCodes/ subdirectory should already exist in the multiPhATE distribution, with PHANOTATE for Python 3.x in that location.)

CONDA INSTALLATION

If you prefer to run multiPhATE in a Conda environment, here are some tips for how to set it up. 

1) First, download and install miniconda3 for Python 3.7 (https://conda.io/en/latest/miniconda.html). For more information about Conda, see https://docs.conda.io/projects/conda/en/latest/user-guide/getting-started.html.  

2) Check that conda is working:  $ conda --version
 
    If conda is not recognized, then you may need to  switch to bash shell: $ bash  
    (and then try again)

3) Add the bioconda channel:  $ conda --add channels bioconda 

    Note: bioconda is supported on Linux and Mac operating systems, but so far not on PC.

4) Create an environment for using multiPhATE; let's call it "multiphate":  $ conda create --name multiphate

5) Activate that environment:  $ source activate multiphate

6) Install conda packages within that environment:  $ conda install python=3

    Repeat for each of biopython, emboss, blast, glimmer, prodigal, hmmer, trnascan-se.  

7) When running multiPhATE within your multiphate Conda environment, the pipeline will use the version of python and the third party codes installed within the multiphate environment, so there should be no clashes with other versions of these packages that may be installed elsewhere on your system. When you are finished running multiPhATE, you may exit from the multiphate Conda environment:  $ source deactivate

Note that genemarks and phanotate are not available as conda packages, so these programs, as well as the dependent databases, all need to be acquired/installed manually in any case.

PHATE PIPELINE OUTPUT FILES 

In the user-specified output directory (<genomeDir>), under PipelineOutput/, the following files will be written:

gene-call outputs from each of the gene callers that was run, including a gff-formatted output file

gene.fnt and protein.faa fasta files generated using the designated preferred gene finder

CGC_results.txt - side-by-side comparison of all gene finder results (if at least two were run) plus comparison statistics

cgc.gff - a superset of gene calls, each tagged with the gene caller(s) that made the call

phate_sequenceAnnotation_main.out - tabbed integrated annotation results (also written in gff format)

BLAST/ directory - raw blast results for genome (under Genome/) and proteins (under Protein/) (if CLEAN_RAW_DATA switch in multiPhate.config file is 'false')

pVOG groupings comprising alignment-ready fasta files (under BLAST/Protein/) based on best blast hits

HMM/ directory - raw hmm search results for protein sequences (under Protein/) (if CLEAN_RAW_DATA switch in phate_runPipeline.py is 'False')

pVOG grouping comprising alignment-ready fasta files (under HMM/Protein/) based on best hmm hits

log files capturing details of the processing and time stamps

The auto-generated myGenomeName_phate.config file, to record exactly how you configured the pipeline for the current run (genome).


RUNNING PHATE AS AN "EMBARASSINGLY PARALLEL" CODE

Pipeline outputs are written to user-specified output subdirectories (specified in your multiPhate.config file, and in the auto-generated myGenomeName_phate.config files, one for each genome). In this way, multiPhATE may be set up to run in parallel on any number of nodes of a compute cluster without fear of clashes in writing results. Implementing a code in parallel is system dependent; a parallel version of multiPhATE is not provided in the distribution. In the case that you have programming expertise at hand, parallelization should be implemented within the multiPhate.py script, upon execution of phate_runPipeline.py over the list of genomes. 


FURTHER RECOMMENDATIONS

PhATE was originally developed for anotating phage genome sequences. However, PhATE may also be useful in helping to identify phage genes within bacterial genomes (i.e., prophage). Thus, the user has the option of running multiple gene callers for bacterial genome sequence (GeneMarkS, Glimmer, Prodigal) and a new gene finder specifically for phage (PHANOTATE). The calls from any two or more of these callers are compared by the PhATE/CGC code so that the user can examine the calls that agree or disagree among the callers, and then run PhATE again, selecting the caller of choice for annotating the sequence.

Although most supporting databases for multiPhATE are phage- or virus-centric, NR and Refseq Protein are included in order to help identify genes/functions that are not known in the virus/phage gene data sets. However, PhATE is not intended to be the sole source for anntotation of bacterial genome sequences, as PhATE is tailored for identification of genes and functions in phage.

The NR database has grown enormously large. It is recommended to use a smaller database, such as Refseq Protein instead of NR. Furthermore, annotating with NR will add greatly to the time required for processing a genome through PhATE. Therefore, it is recommended that NR be turned off ('false') until one desires to preform a full/final annotation of the genome of interest, if using NR.

Because the behavior of 3rd party codes can sometimes be unpredictable, it is recommended that the user replace spaces and unusual characters in their fasta headers with the underscore character.


PLANNED FURTHER DEVELOPMENT

Plans include adding additional functionality to the pipeline:
1) Adding HMMER and additional hmm databases
2) Implementing an optional "custom" database for blasting or doing hmm analysis. (Actually, you can do custom blast now by substituting your custom nucleotide blast database for the refseq_gene_blast process, or substituting your custom protein blast database for one of the protein processes (swissprot, nr, refseq_protein, or ncbi_virus_protein), by specifying the location of your custom database in the multiPhate.config file and turning the corresponding annotation process to 'true'.)

Feel free to report bugs or problems, or to suggest future improvements, by posting an issue on the project github page (click on the Issues tab), or by emailing the developers at:  zhou4@llnl.gov. Thank you for using multiPhATE.

PUBLICATION

If you use multiPhATE in your research, kindly reference our paper:  "multiPhATE: bioinformatics pipeline for functional annotation of phage isolates", by Carol E Zhou, Stephanie A Malfatti, Jeffrey A Kimbrel, Casandra W Philipson, Katelyn E McNair, Theron C Hamilton, Robert A Edwards, and Brian E Souza. BioRxiv Feb 15, 2019, doi: http://dx.doi.org/10.1101/551010. 

multiPhATE v.0.5
