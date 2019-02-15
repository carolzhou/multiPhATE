#!/usr/bin/env python

################################################################
#
# Program Title:  multiPhate.py (/MultiPhate/)
#
# Last Update:  14 February 2019
#
# Description: Script multiPhate.py runs the phate annotation pipeline over a set of input phage genomes.  This code runs under 
#    Python 2.7, and requires dependent packages and databases as listed in the README file.
#    multiPhate.py inputs a configuration file (see sample_ multiPhate.config), and uses it to construct a set of
#    configuration files, one for each genome. Then, multiPhate.py executes phate_runPipeline.py over all of the genomes in the set.
#
# Usage:  python multiPhate.py myMultiPhate.config
#    (see sample_multiPhate.config for how to create your configuration file)
#
################################################################

# This code was developed by Carol L. Ecale Zhou at Lawrence Livermore National Laboratory.
# THIS CODE IS COVERED BY THE BSD LICENSE. SEE INCLUDED FILE BSD.PDF FOR DETAILS.

# DO NOT MODIFY ANYTHING IN THIS FILE EXCEPT ITEMS LABELED AS "USER CONFIGURATION"

import sys, os, re, string, copy, time, datetime
import subprocess

# CONFIGURABLE
# 1) If you are running under a linux system, set PHATE_OUT and PHATE_ERR to 'True'. This will capture standard errors to files. Cannot
# guarantee this will work under other operating systems.
PHATE_OUT = 'False'
PHATE_ERR = 'True'
#
# Default Verbosity
CLEAN_RAW_DATA_DEFAULT = 'True'   # if 'False', the raw Blast and Hmm outputs will be saved in the PipelineOutput folder
PHATE_WARNINGS_DEFAULT = 'False'
PHATE_MESSAGES_DEFAULT = 'False'
PHATE_PROGRESS_DEFAULT = 'True'
CGC_WARNINGS_DEFAULT   = 'False'
CGC_MESSAGES_DEFAULT   = 'False'
CGC_PROGRESS_DEFAULT   = 'False'
#DEBUG = True     # Controls debug settings in this (local) code only
DEBUG = False    # Leave False, unless debugging

# Env:  BLAST parameters. Normally leave these settings alone. These are minimum cutoffs. Configure stringency in config file.
BLASTP_IDENTITY_DEFAULT  = '60'
BLASTP_HIT_COUNT_DEFAULT = '3'
BLASTN_HIT_COUNT_DEFAULT = '3'

# Constants; defaults will apply if not specified in config file
# Leave all this stuff alone! 

# Standard directories
BASE_DIR_DEFAULT            = os.path.join(os.getcwd(),"")    # Ex: /Home/MyName/MyCodeDirectory/multiPhATE/
DATABASE_DIR_DEFAULT        = BASE_DIR_DEFAULT + "Databases/" 
SOFTWARE_DIR_DEFAULT        = BASE_DIR_DEFAULT + "ExternalCodes/" 
PIPELINE_INPUT_DIR_DEFAULT  = BASE_DIR_DEFAULT + "PipelineInput/"    
PIPELINE_OUTPUT_DIR_DEFAULT = BASE_DIR_DEFAULT + "PipelineOutput/"   

PHATE_PIPELINE_CODE         = 'phate_runPipeline.py'
CONSENSUS_CALLS_FILE        = 'phanotate.cgc' #*** For now this is PHANOTATE calls, though may be consensus calls in future
GENE_FILE                   = 'gene.fnt'      #
PROTEIN_FILE                = 'protein.faa'   #
GENETIC_CODE                = '11'        # default is bacterial (11)
GENE_CALLER                 = 'phanotate' # default is annotation of phage, so PHANOTATE is preferred gene caller; if bac, could be 'consensus', 'genemark', 'glimmer', or 'prodigal'
GENOME_TYPE                 = 'phage'     # default is phage; could be 'bacterium'
NAME                        = 'unknown'   # user provided
CONTIG_NAME                 = 'unknown'   # user provided: temporary, finished genomes/single contig only for now
SPECIES                     = 'unknown'   # user provided
 
# gene callers
GENEMARKS_CALLS_DEFAULT     = False     # Requires license
PRODIGAL_CALLS_DEFAULT      = False
GLIMMER_CALLS_DEFAULT       = False
PHANOTATE_CALLS_DEFAULT     = False

GENEMARKS_PATH              = ''        # Available via license	 
GLIMMER_PATH                = ''        # Can install using Conda
PRODIGAL_PATH               = ''        # Can install using Conda
PHANOTATE_PATH              = ''

#blast parameters
MAX_BLAST_HIT_COUNT         = 100     # maximum number of hits to capture (user should specify far fewer than max)
MIN_BLASTP_IDENTITY         = 20      # default; sets a lower limit based on value at which a structure model can provide information
MAX_BLASTP_HIT_COUNT        = 100     # default; sets an upper limit; user's value should typically be well below this
MAX_BLASTN_HIT_COUNT        = 10      # default; sets an upper limit

#blast databases to be used for search
NCBI_VIRUS_BLAST_DEFAULT         = False
NCBI_VIRUS_PROTEIN_BLAST_DEFAULT = False
KEGG_VIRUS_BLAST_DEFAULT         = False     # Requires license
NR_BLAST_DEFAULT                 = False     # Large data set; blast run takes time
REFSEQ_PROTEIN_BLAST_DEFAULT     = False     # Large data set; blast run takes time
PHANTOME_BLAST_DEFAULT           = False
PVOGS_BLAST_DEFAULT              = False
UNIPARC_BLAST_DEFAULT            = False     # Keep turned 'off' for now; not yet in service
REFSEQ_GENE_BLAST_DEFAULT        = False
SWISSPROT_BLAST_DEFAULT          = False
UNIPROT_BLAST_DEFAULT            = False     # not yet in service
PFAM_BLAST_DEAFULT               = False     # not yet in service

#hmm programs
HMM_PROGRAM_DEFAULT              = 'jackhmmer'  # This is the only hmm program currently supported

#hmm databases to be used for search
NCBI_VIRUS_HMM_DEFAULT           = False     # not yet in service
NCBI_VIRUS_PROTEIN_HMM_DEFAULT   = False     # not yet in service
KEGG_VIRUS_HMM_DEFAULT           = False     # Requires license
NR_HMM_DEFAULT                   = False     # Large data set; hmm run takes time
REFSEQ_PROTEIN_HMM_DEFAULT       = False     # Large data set; hmm run takes time
PHANTOME_HMM_DEFAULT             = False     # not yet in service
PVOGS_HMM_DEFAULT                = False     #
UNIPARC_HMM_DEFAULT              = False     # not yet in service
REFSEQ_GENE_HMM_DEFAULT          = False     # not yet in service
SWISSPROT_HMM_DEFAULT            = False     # not yet in service
UNIPROT_HMM_DEFAULT              = False     # not yet in service
PFAM_HMM_DEAFULT                 = False     # not yet in service

#other
PSAT_ANNOTATION_DEFAULT          = False     # Requires LLNL processing
PSAT                             = False
PSAT_FILE                        = ""

# ENVIRONMENT VARIABLES
# It is most convenient to locate the supporting software codes and databases in the above-indicated subdirectories.
# However, if any of your supporting databases or softwares reside elsewhere, then explicit locations will need to 
# be filled in in the multiPhate.config file. This will likely be the case for large databases that you may already 
# have on your compute cluster (e.g, NR), and for software packages, such as EMBOSS or gene finders that you may
# already have installed on your system. Parameters that differ from defaults will be re-assigned based on information
# provided in the users' multiPhate.config file.

PIPELINE_INPUT_DIR                          = BASE_DIR_DEFAULT + PIPELINE_INPUT_DIR_DEFAULT   # Default
PIPELINE_OUTPUT_DIR                         = BASE_DIR_DEFAULT + PIPELINE_OUTPUT_DIR_DEFAULT  # Default
PHATE_BASE_DIR                              = BASE_DIR_DEFAULT
EMBOSS_CODE                                 = ""  # Modify this for the version you have
EMBOSS_PHATE_HOME                           = SOFTWARE_DIR_DEFAULT + EMBOSS_CODE     # if installed in SOFTWARE_DIR, else enter actual location
os.environ["BASE_DIR"]                      = BASE_DIR_DEFAULT  
os.environ["DATABASE_DIR"]                  = DATABASE_DIR_DEFAULT
os.environ["SOFTWARE_DIR"]                  = SOFTWARE_DIR_DEFAULT
os.environ["PIPELINE_INPUT_DIR"]            = PIPELINE_INPUT_DIR
os.environ["PIPELINE_OUTPUT_DIR"]           = PIPELINE_OUTPUT_DIR
os.environ["PHATE_BASE_DIR"]                = PHATE_BASE_DIR
os.environ["EMBOSS_PHATE_HOME"]             = EMBOSS_PHATE_HOME 
os.environ["PIPELINE_DIR"]                  = BASE_DIR_DEFAULT
os.environ["PSAT_OUT_DIR"]                  = BASE_DIR_DEFAULT

# Data sets
os.environ["KEGG_VIRUS_BASE_DIR"]           = DATABASE_DIR_DEFAULT + "KEGG/"
os.environ["KEGG_VIRUS_BLAST_HOME"]         = os.environ["KEGG_VIRUS_BASE_DIR"] + "T40000.pep"
os.environ["NCBI_VIRUS_BASE_DIR"]           = DATABASE_DIR_DEFAULT + "NCBI/"
os.environ["NCBI_VIRUS_BLAST_HOME"]         = os.environ["NCBI_VIRUS_BASE_DIR"] + "Virus_Genome/"  + "viral.1.1.genomic.fna"
os.environ["NCBI_VIRUS_PROTEIN_BLAST_HOME"] = os.environ["NCBI_VIRUS_BASE_DIR"] + "Virus_Protein/" + "viral.protein.faa"
os.environ["NCBI_TAXON_DIR"]                = os.environ["NCBI_VIRUS_BASE_DIR"] + "Virus_Genome/"
os.environ["PHANTOME_BASE_DIR"]             = DATABASE_DIR_DEFAULT + "Phantome/"
os.environ["PHANTOME_BLAST_HOME"]           = os.environ["PHANTOME_BASE_DIR"] + "Phantome_Phage_genes.faa"
os.environ["PVOGS_BASE_DIR"]                = DATABASE_DIR_DEFAULT + "pVOGs/"
os.environ["PVOGS_BLAST_HOME"]              = os.environ["PVOGS_BASE_DIR"] + "pVOGs.faa"
os.environ["UNIPARC_BASE_DIR"]              = DATABASE_DIR_DEFAULT + "UniParc/" # Uniparc not yet in service
os.environ["UNIPARC_VIRUS_BLAST_HOME"]      = os.environ["UNIPARC_BASE_DIR"] + "uniparc_active.fasta"  #*** ???
os.environ["NR_BLAST_BASE_DIR"]             = "/data/data1/sandbox/BLAST/"
os.environ["NR_BLAST_HOME"]                 = os.environ["NR_BLAST_BASE_DIR"] + "nr"
os.environ["REFSEQ_PROTEIN_BASE_DIR"]       = DATABASE_DIR_DEFAULT + "Refseq/Protein/"
os.environ["REFSEQ_PROTEIN_BLAST_HOME"]     = os.environ["REFSEQ_PROTEIN_BASE_DIR"] + "refseq_protein"
os.environ["REFSEQ_GENE_BASE_DIR"]          = DATABASE_DIR_DEFAULT + "Refseq/Gene/"
os.environ["REFSEQ_GENE_BLAST_HOME"]        = os.environ["REFSEQ_GENE_BASE_DIR"] + "refseqgene"
os.environ["SWISSPROT_BASE_DIR"]            = DATABASE_DIR_DEFAULT + "Swissprot/"
os.environ["SWISSPROT_BLAST_HOME"]          = os.environ["SWISSPROT_BASE_DIR"] + "swissprot"
os.environ["UNIPROT_BASE_DIR"]              = DATABASE_DIR_DEFAULT + "Uniprot/" # not yet in service
os.environ["UNIPROT_BLAST_HOME"]            = os.environ["UNIPROT_BASE_DIR"] + "uniprot"
os.environ["PFAM_BASE_DIR"]                 = DATABASE_DIR_DEFAULT + "Pfam/" # not yet in service
os.environ["PFAM_BLAST_HOME"]               = os.environ["PFAM_BASE_DIR"] + "pfam"

# Gene calling
#os.environ["PRODIGAL_PATH"]                 = SOFTWARE_DIR_DEFAULT + "prodigal.v2_50/"
os.environ["PRODIGAL_PATH"]                 = ""   # global, if installed via conda
#os.environ["GLIMMER_PATH"]                  = SOFTWARE_DIR_DEFAULT + "glimmer3.02/bin/"
os.environ["GLIMMER_PATH"]                  = ""   # global, if installed via conda
os.environ["GENEMARKS_PATH"]                = SOFTWARE_DIR_DEFAULT + "GeneMarkS/genemark_suite_linux_64/gmsuite/"
os.environ["PHANOTATE_PATH"]                = SOFTWARE_DIR_DEFAULT + "PHANOTATE/PHANOTATE-master/"
#os.environ["PHANOTATE_PATH"]                = ""   # global, if installed via conda
os.environ["CGC_PATH"]                      = BASE_DIR_DEFAULT + "CompareCalls/"
#os.environ["tRNAscanSE_HOME"]               = /Users/myName/tRNAscanDir/trnascan-se-2.0/bin/tRNAscan-SE"
os.environ["tRNAscanSE_HOME"]               = ""    # global, if installed via conda

# Blast
#os.environ["BLAST_HOME"]                    = SOFTWARE_DIR_DEFAULT + "/ncbi-blast-2.7.1+/bin/"
os.environ["BLAST_HOME"]                    = ""    # global, if installed via conda; use blast+, not legacy blast 
os.environ["MIN_BLASTP_IDENTITY"]           = str(MIN_BLASTP_IDENTITY)
os.environ["MAX_BLASTP_HIT_COUNT"]          = str(MAX_BLASTP_HIT_COUNT)
os.environ["MAX_BLASTN_HIT_COUNT"]          = str(MAX_BLASTN_HIT_COUNT)
os.environ["BLASTP_IDENTITY_DEFAULT"]       = str(BLASTP_IDENTITY_DEFAULT)
os.environ["BLASTP_HIT_COUNT_DEFAULT"]      = str(BLASTP_HIT_COUNT_DEFAULT)
os.environ["BLASTN_HIT_COUNT_DEFAULT"]      = str(BLASTN_HIT_COUNT_DEFAULT)

# HMM
os.environ["HMM_HOME"]                      = ""

# Global control: verbosity and error capture
os.environ["CLEAN_RAW_DATA"]                = CLEAN_RAW_DATA_DEFAULT
os.environ["PHATE_WARNINGS"]                = PHATE_WARNINGS_DEFAULT  # Print warnings and errors to standard out
os.environ["PHATE_MESSAGES"]                = PHATE_MESSAGES_DEFAULT  # Print helpful messages (may be verbose)
os.environ["PHATE_PROGRESS"]                = PHATE_PROGRESS_DEFAULT  # Print each step in processing a genome
os.environ["PHATE_ERR"]                     = PHATE_ERR       # Capture standard errors to files on linux/mac machine
os.environ["PHATE_OUT"]                     = PHATE_OUT       # Capture standard errors to files on linux/mac machine
os.environ["CGC_WARNINGS"]                  = CGC_WARNINGS_DEFAULT
os.environ["CGC_MESSAGES"]                  = CGC_MESSAGES_DEFAULT
os.environ["CGC_PROGRESS"]                  = CGC_PROGRESS_DEFAULT


# Constants

CODE_BASE   = "multiPhate"
CODE        = CODE_BASE + ".py"
CONFIG_FILE = "multiPhate.config"            # by default, but user should name their own, ending in ".config"
SAMPLE_CONFIG_FILE = "sample_" + CONFIG_FILE # Sample config file; user should copy, then modify. 

# HELP STRINGS

HELP_STRING = """This code, """ + CODE + """, runs the phage annotation pipeline (phate_runPipeline.py) over multipe genomes. The configuration file input to this code specifies a list of genomes to be processed and the parameters for pipeline execution. The pipeline performs 1) gene calling by 4 gene callers (PHANOTATE, GeneMarkS, Glimmer3, and Prodigal), followed by identification of closest phage genome by means of blast against an NCBI-phage database, and sequence-based functional annotation by means of blastp against several peptide databases (NR, NCBI virus protein, KEGG-virus, Phantome, pVOGs, Swissprot, Refseq protein), and HMM search against the pVOG database. \nType: python """ + CODE + """ usage - for more information about constructing the command line.\nType: python """ + CODE + """ detail - for more information about how this code can be run.\n"""

INPUT_STRING = """The input files and other parameters for running this code are specified in a configuration file, which is provided as the only input parameter. See sample configuration file (""" + SAMPLE_CONFIG_FILE + """) for details on how to customize your configuration file.\n"""

USAGE_STRING = """Usage: python """ + CODE + """ """ + CONFIG_FILE + """\n"""

DETAIL_STRING = """Currently the PSAT module is run separately as a web service. In order to incorporate PSAT output into your annotations, you should first run this pipeline specifying "translation_only" in the configuration file. Then, use the generated peptide/protein fasta file as input for PSAT processing. Once you have the PSAT output, save it to the pipeline input directory, and re-run this pipeline, specifying that translation_only is false.\n"""

##### PATTERNS #####

# Locations
p_phateDir                    = re.compile("phate_dir='(.*)'")
p_databaseDir                 = re.compile("database_dir='(.*)'")
p_softwareDir                 = re.compile("software_dir='(.*)'")

# General
p_comment                     = re.compile("^#")
p_blank                       = re.compile("^$")
p_help                        = re.compile("help")
p_input                       = re.compile("input")
p_usage                       = re.compile("usage")
p_detail                      = re.compile("detail")
p_config                      = re.compile("config")
p_outputSubdir                = re.compile("output_subdir='(.*)'")
p_genomeFile                  = re.compile("genome_file='(.*)'")
p_genomeType                  = re.compile("genome_type='(.*)'")
p_name                        = re.compile("name='(.*)'")
p_contig                      = re.compile("contig='(.*)'")  #*** For now, finished genome, single contig only
p_species                     = re.compile("species='(.*)'")

# Genome information
p_genomeList                  = re.compile("Genome\sList")  # non-case-sensitive "Genome List"
p_genomeNumber                = re.compile("Genome\s+(\d+)")     # genome number
p_root                        = re.compile("([\w\d_-]+)\.fasta")    # captures the root name of the fasta file (e.g., takes 'P2' from P2.fasta)
p_end                         = re.compile("END")

# Gene calling
p_geneCaller                  = re.compile("gene_caller='(.*)'")
p_genemarksCalls              = re.compile("genemarks_calls='(.*)'")
p_glimmerCalls                = re.compile("glimmer_calls='(.*)'")
p_prodigalCalls               = re.compile("prodigal_calls='(.*)'")
p_phanotateCalls              = re.compile("phanotate_calls='(.*)'")
p_geneticCode                 = re.compile("genetic_code='(\d+)'")
p_translateOnly               = re.compile("translate_only='(.*)'")

# Blast
p_blastpIdentity              = re.compile("blast_identity='(\d+)'")   #*** For now; but should distinguish between blastn/blastp
p_blastpHitCount              = re.compile("blastp_hit_count='(\d+)'")
p_blastnHitCount              = re.compile("blastn_hit_count='(\d+)'")
p_ncbiVirusBlast              = re.compile("ncbi_virus_blast='(.*)'")
p_ncbiVirusProteinBlast       = re.compile("ncbi_virus_protein_blast='(.*)'")
p_keggVirusBlast              = re.compile("kegg_virus_blast='(.*)'")
p_nrBlast                     = re.compile("nr_blast='(.*)'")
p_refseqProteinBlast          = re.compile("refseq_protein_blast='(.*)'")
p_refseqGeneBlast             = re.compile("refseq_gene_blast='(.*)'")
p_swissprotBlast              = re.compile("swissprot_blast='(.*)'")
p_phantomeBlast               = re.compile("phantome_blast='(.*)'")
p_pvogsBlast                  = re.compile("pvogs_blast='(.*)'")
p_uniparcBlast                = re.compile("uniparc_blast='(.*)'")
p_uniprotBlast                = re.compile("uniprot_blast='(.*)'")
p_refseqGeneBlast             = re.compile("refseq_gene_blast='(.*)'")

# HMM
p_hmmProgram                  = re.compile("hmm_program='(.*)'")
p_ncbiVirusHmm                = re.compile("ncbi_virus_hmm='(.*)'")
p_ncbiVirusProteinHmm         = re.compile("ncbi_virus_protein_hmm='(.*)'")
p_keggVirusHmm                = re.compile("kegg_virus_hmm='(.*)'")
p_phantomeHmm                 = re.compile("phantome_hmm='(.*)'")
p_pvogsHmm                    = re.compile("pvogs_hmm='(.*)'")
p_swissprotHmm                = re.compile("swissprot_hmm='(.*)'")
p_refseqProteinHmm            = re.compile("refseq_protein_hmm='(.*)'")
p_refseqGeneHmm               = re.compile("refseq_gene_hmm='(.*)'")
p_uniparcHmm                  = re.compile("uniparc_hmm='(.*)'")
p_uniprotHmm                  = re.compile("uniprot_hmm='(.*)'")
p_nrHmm                       = re.compile("nr_hmm='(.*)'")

# Dependent Code Locations
p_blastPlusHome               = re.compile("blast_plus_home='(.*)'")
p_embossHome                  = re.compile("emboss_home='(.*)'")
p_tRNAscanSEhome              = re.compile("tRNAscanSE_home='(.*)'")
p_glimmerHome                 = re.compile("glimmer_home='(.*)'")
p_prodigalHome                = re.compile("prodigal_home='(.*)'")
p_phanotateHome               = re.compile("phanotate_home='(.*)'")
p_genemarkHome                = re.compile("genemark_home='(.*)'")

# Database Locations
p_ncbiVirusDatabase           = re.compile("ncbi_virus_database='(.*)'")
p_refseqGeneDatabase          = re.compile("refseq_gene_database='(.*)'")
p_ncbiVirusProteinDatabase    = re.compile("ncbi_virus_protein_database='(.*)'")
p_keggVirusDatabase           = re.compile("kegg_virus_database='(.*)'")
p_phantomeDatabase            = re.compile("phantome_database='(.*)'")
p_pvogsDatabase               = re.compile("pvogs_database='(.*)'")
p_swissprotDatabase           = re.compile("swissprot_database='(.*)'")
p_refseqProteinDatabase       = re.compile("refseq_protein_database='(.*)'")
p_nrDatabase                  = re.compile("nr_database='(.*)'")

# Verbosity
p_phateWarnings               = re.compile("phate_warnings='(.*)'")
p_phateMessages               = re.compile("phate_messages='(.*)'")
p_phateProgress               = re.compile("phate_progress='(.*)'")
p_cgcWarnings                 = re.compile("cgc_warnings='(.*)'")
p_cgcMessages                 = re.compile("cgc_messages='(.*)'")
p_cgcProgress                 = re.compile("cgc_progress='(.*)'")
p_cleanRawData                = re.compile("clean_raw_data='(.*)'")

# PSAT
p_psatAnnotation              = re.compile("psat_annotation='(.*)'")
p_psatFile                    = re.compile("psat_file='(.*)'")

##### GET INPUT PARAMETERS #####

# Open log file
logfile = BASE_DIR_DEFAULT + CODE_BASE + ".log"
LOG = open(logfile,'w')
LOG.write("%s%s\n" % ("Begin log file ",datetime.datetime.now()))

if len(sys.argv) != 2:
    print(HELP_STRING)
    dateTime = os.popen('date')
    LOG.write("%s%s%s%s\n" % ("Incorrect number of input parameters: ", len(sys.argv), ". End log ",dateTime))
    LOG.close(); exit(0)
else:
    match_config = re.search(p_config,sys.argv[1])
    if match_config:
        configFile = sys.argv[1]
        LOG.write("%s%s\n" % ("Config file is ",configFile))
    else: 
        match_input  = re.search(p_input,  sys.argv[1].lower())
        match_usage  = re.search(p_usage,  sys.argv[1].lower())
        match_detail = re.search(p_detail, sys.argv[1].lower())
        if match_input:
            print(INPUT_STRING)
        elif match_usage:
            print(USAGE_STRING)
        elif match_detail:
            print(DETAIL_STRING)
        else:
            print(HELP_STRING)
        LOG.write("%s%s\n" % ("A help string was provided to user; End log ",datetime.datetime.now()))
        LOG.close(); exit(0)

# Open and check input file

fileError = False
try:
    CONFIG = open(configFile,"r")
except IOError as e:
    fileError = True
    print(e)

if fileError:
    print("Check your config file.")
    print(HELP_STRING)
    LOG.write("%s%s\n" % ("A help string was provided to user; End log ",datetime.datetime.now()))
    LOG.close(); exit(0)

##### Read input parameters from configuration file

# First, set as defaults; note: setting these values in config file is optional

geneticCode           = GENETIC_CODE
geneCaller            = GENE_CALLER
genomeType            = GENOME_TYPE
name                  = NAME
contigName            = CONTIG_NAME
species               = SPECIES

blastpIdentity        = BLASTP_IDENTITY_DEFAULT
blastpHitCount        = BLASTP_HIT_COUNT_DEFAULT
blastnHitCount        = BLASTN_HIT_COUNT_DEFAULT
ncbiVirusBlast        = NCBI_VIRUS_BLAST_DEFAULT
ncbiVirusProteinBlast = NCBI_VIRUS_PROTEIN_BLAST_DEFAULT
keggVirusBlast        = KEGG_VIRUS_BLAST_DEFAULT
nrBlast               = NR_BLAST_DEFAULT
refseqProteinBlast    = REFSEQ_PROTEIN_BLAST_DEFAULT
refseqGeneBlast       = REFSEQ_GENE_BLAST_DEFAULT
phantomeBlast         = PHANTOME_BLAST_DEFAULT
pvogsBlast            = PVOGS_BLAST_DEFAULT
uniparcBlast          = UNIPARC_BLAST_DEFAULT
uniprotBlast          = UNIPROT_BLAST_DEFAULT
swissprotBlast        = SWISSPROT_BLAST_DEFAULT

hmmProgram            = HMM_PROGRAM_DEFAULT
ncbiVirusHmm          = NCBI_VIRUS_HMM_DEFAULT
ncbiVirusProteinHmm   = NCBI_VIRUS_PROTEIN_HMM_DEFAULT
keggVirusHmm          = KEGG_VIRUS_HMM_DEFAULT
nrHmm                 = NR_HMM_DEFAULT
refseqGeneHmm         = REFSEQ_GENE_HMM_DEFAULT
refseqProteinHmm      = REFSEQ_PROTEIN_HMM_DEFAULT
phantomeHmm           = PHANTOME_HMM_DEFAULT
pvogsHmm              = PVOGS_HMM_DEFAULT
uniparcHmm            = UNIPARC_HMM_DEFAULT
uniprotHmm            = UNIPROT_HMM_DEFAULT
swissprotHmm          = SWISSPROT_HMM_DEFAULT

genemarksCalls        = GENEMARKS_CALLS_DEFAULT
prodigalCalls         = PRODIGAL_CALLS_DEFAULT
glimmerCalls          = GLIMMER_CALLS_DEFAULT
phanotateCalls        = PHANOTATE_CALLS_DEFAULT
psatAnnotation        = PSAT_ANNOTATION_DEFAULT

# Capture user's configured values

FIRST_GENOME = True  
DATA_ITEMS_NUM = 6 
genomeDataDict = {
    "genomeNumber"  : "",
    "genomeFile"    : "",
    "genomeType"    : "",
    "genomeSpecies" : "",
    "genomeName"    : "",
    "outputSubdir"  : "",
    }
genomeList      = []  # List of genomeData objects
genomeNumber    = ""  # Number of current genome; assigned by user; could be a string
genomeDataItems = 0   # Number of data items collected for current genome; should be DATA_ITEMS_NUM 
BEGIN_GENOME_LIST = False
nextGenomeData = genomeDataDict

cLines = CONFIG.read().splitlines()
for cLine in cLines:
    match_comment                  = re.search(p_comment,cLine)
    match_blank                    = re.search(p_blank,cLine)

    match_genomeList               = re.search(p_genomeList,cLine)
    match_genomeNumber             = re.search(p_genomeNumber,cLine)
    match_genomeFile               = re.search(p_genomeFile,cLine)
    match_genomeType               = re.search(p_genomeType,cLine)
    match_species                  = re.search(p_species,cLine)
    match_name                     = re.search(p_name,cLine)
    match_outputSubdir             = re.search(p_outputSubdir,cLine)
    match_end                      = re.search(p_end,cLine)

    match_psatFile                 = re.search(p_psatFile,cLine)
    match_geneticCode              = re.search(p_geneticCode,cLine)
    match_translateOnly            = re.search(p_translateOnly,cLine)
    match_geneCaller               = re.search(p_geneCaller,cLine)
    match_contig                   = re.search(p_contig,cLine)

    match_blastpIdentity           = re.search(p_blastpIdentity,cLine)
    match_blastpHitCount           = re.search(p_blastpHitCount,cLine)
    match_blastnHitCount           = re.search(p_blastnHitCount,cLine)
    match_ncbiVirusBlast           = re.search(p_ncbiVirusBlast,cLine)
    match_ncbiVirusProteinBlast    = re.search(p_ncbiVirusProteinBlast,cLine)
    match_keggVirusBlast           = re.search(p_keggVirusBlast,cLine)
    match_nrBlast                  = re.search(p_nrBlast,cLine)
    match_refseqProteinBlast       = re.search(p_refseqProteinBlast,cLine)
    match_refseqGeneBlast          = re.search(p_refseqGeneBlast,cLine)
    match_phantomeBlast            = re.search(p_phantomeBlast,cLine)
    match_pvogsBlast               = re.search(p_pvogsBlast,cLine)
    match_uniparcBlast             = re.search(p_uniparcBlast,cLine)
    match_uniprotBlast             = re.search(p_uniprotBlast,cLine)
    match_swissprotBlast           = re.search(p_swissprotBlast,cLine)
    match_refseqGeneBlast          = re.search(p_refseqGeneBlast,cLine)

    match_hmmProgram               = re.search(p_hmmProgram,cLine)
    match_ncbiVirusHmm             = re.search(p_ncbiVirusHmm,cLine)
    match_ncbiVirusProteinHmm      = re.search(p_ncbiVirusProteinHmm,cLine)
    match_keggVirusHmm             = re.search(p_keggVirusHmm,cLine)
    match_nrHmm                    = re.search(p_nrHmm,cLine)
    match_refseqProteinHmm         = re.search(p_refseqProteinHmm,cLine)
    match_phantomeHmm              = re.search(p_phantomeHmm,cLine)
    match_pvogsHmm                 = re.search(p_pvogsHmm,cLine)
    match_uniparcHmm               = re.search(p_uniparcHmm,cLine)
    match_uniprotHmm               = re.search(p_uniprotHmm,cLine)
    match_swissprotHmm             = re.search(p_swissprotHmm,cLine)
    match_refseqGeneHmm            = re.search(p_refseqGeneHmm,cLine)

    match_genemarksCalls           = re.search(p_genemarksCalls,cLine)
    match_prodigalCalls            = re.search(p_prodigalCalls,cLine)
    match_glimmerCalls             = re.search(p_glimmerCalls,cLine)
    match_phanotateCalls           = re.search(p_phanotateCalls,cLine)
    match_psatAnnotation           = re.search(p_psatAnnotation,cLine)

    match_phateDir                 = re.search(p_phateDir,cLine)
    match_databaseDir              = re.search(p_databaseDir,cLine)
    match_softwareDir              = re.search(p_softwareDir,cLine)

    match_blastPlusHome            = re.search(p_blastPlusHome,cLine)
    match_embossHome               = re.search(p_embossHome,cLine)
    match_tRNAscanSEhome           = re.search(p_tRNAscanSEhome,cLine)
    match_glimmerHome              = re.search(p_glimmerHome,cLine)
    match_prodigalHome             = re.search(p_prodigalHome,cLine)
    match_phanotateHome            = re.search(p_phanotateHome,cLine)
    match_genemarkHome             = re.search(p_genemarkHome,cLine)

    match_ncbiVirusDatabase        = re.search(p_ncbiVirusDatabase,cLine)
    match_refseqGeneDatabase       = re.search(p_refseqGeneDatabase,cLine)
    match_ncbiVirusProteinDatabase = re.search(p_ncbiVirusProteinDatabase,cLine)
    match_keggVirusDatabase        = re.search(p_keggVirusDatabase,cLine)
    match_phantomeDatabase         = re.search(p_phantomeDatabase,cLine)
    match_pvogsDatabase            = re.search(p_pvogsDatabase,cLine)
    match_swissprotDatabase        = re.search(p_swissprotDatabase,cLine)
    match_refseqProteinDatabase    = re.search(p_refseqProteinDatabase,cLine)
    match_nrDatabase               = re.search(p_nrDatabase,cLine)

    match_phateWarnings            = re.search(p_phateWarnings,cLine)
    match_phateMessages            = re.search(p_phateMessages,cLine)
    match_phateProgress            = re.search(p_phateProgress,cLine)
    match_cgcWarnings              = re.search(p_cgcWarnings,cLine)
    match_cgcMessages              = re.search(p_cgcMessages,cLine)
    match_cgcProgress              = re.search(p_cgcProgress,cLine)
    match_cleanRawData             = re.search(p_cleanRawData,cLine)
 
    ##### Capture list of genomes and associated data #####

    if (match_comment or match_blank):
        pass 

    elif match_genomeList:  # Capture all genomes listed; for each, gather genome file, genome type, species, name, output subdir
        pass 

    elif match_genomeNumber:  # The next genome's data 
        # First, record previous genome data, if this is not the first genome 
        if not FIRST_GENOME:
            #LOG.write("%s\n" % ("Appending a genome data set"))
            if genomeDataItems != DATA_ITEMS_NUM:  # If record appears incomplete, flag a problem
                LOG.write("%s%s\n" % ("WARNING: check config file for possible incorrect data items: ", genomeDataItems))
            genomeList.append(nextGenomeData)

        # Next, begin collecting next genome's data
        #LOG.write("%s%s\n" % ("Creating a new genome data set for ", match_genomeNumber.group(0)))
        genomeDataItems = 0
        nextGenomeData = copy.deepcopy(genomeDataDict)  # make new genome data object
        genomeNumber = match_genomeNumber.group(1)
        nextGenomeData["genomeNumber"] = genomeNumber
        genomeDataItems += 1
        FIRST_GENOME = False

    elif match_genomeFile:
        value = match_genomeFile.group(1)
        if value != '':
            GENOME_FILE = value 
        else:
            GENOME_FILE = "unknown"
            if PHATE_WARNINGS:
                print("multiPhate says, WARNING:  GENOME_FILE is", GENOME_FILE)
        nextGenomeData["genomeFile"] = GENOME_FILE
        #LOG.write("%s%s\n" % ("GENOME_FILE is ",GENOME_FILE))
        genomeDataItems += 1

    elif match_genomeType:
        value = match_genomeType.group(1)
        if value.lower() == 'phage' or value.lower() == 'bacteriophage':
            nextGenomeData["genomeType"] = 'phage' 
        elif value.lower() == 'virus' or value.lower() == 'viral' or value.lower() == 'viridae':
            nextGenomeData["genomeType"] = 'virus'
        elif value.lower() == 'bacteria' or value.lower() == 'bacterium' or value.lower() == 'bacterial':
            nextGenomeData["genomeType"] = 'bacterium' 
        else:
            nextGenomeData["genomeType"] = 'other'
        #LOG.write("%s%s\n" % ("genome type is ",nextGenomeData["genomeType"]))
        genomeDataItems += 1

    elif match_species:
        species = match_species.group(1)
        nextGenomeData["genomeSpecies"] = species 
        #LOG.write("%s%s\n" % ("Species is ",species))
        genomeDataItems += 1

    elif match_name:
        name = match_name.group(1)
        nextGenomeData["genomeName"] = name 
        #LOG.write("%s%s\n" % ("genome name is ",name))
        genomeDataItems += 1

    elif match_phateDir:
        if match_phateDir.group(1) != '':
            os.environ["BASE_DIR"] = match_phateDir.group(1)

    elif match_databaseDir:
        if match_databaseDir.group(1) != '':
            os.environ["DATABASE_DIR"] = match_databaseDir.group(1)

    elif match_softwareDir:
        if match_softwareDir.group(1) != '':
            os.environ["SOFTWARE_DIR"] = match_softwareDir.group(1)

    elif match_outputSubdir: #*** Note that if the output dir is not read before subdir; depends on user not changing order in config - Clean this up!
        value = match_outputSubdir.group(1)
        if value != '':
            value = value.rstrip('/')  # be sure that name of subdir ends in exactly one '/' (user might omit the slash)
            value = value + '/'
            nextGenomeData["outputSubdir"] = value
        else:
            nextGenomeData["outputSubdir"] = "unknown" 
            if PHATE_WARNINGS:
                print("multiPhate says, WARNING: pipeline output subdir is ", "unknown") 
        genomeDataItems += 1

    elif match_end:  # List of genomes complete; record last genome's data
        if genomeDataItems != DATA_ITEMS_NUM:
            LOG.write("%s%s%s%s%s%s\n" % ("multiPhate says, WARNING: check config file for possible incorrect data items: ", genomeDataItems, " for genome ",nextGenomeData["genomeName"],' ',nextGenomeData["genomeNumber"]))
        genomeList.append(nextGenomeData)
        LOG.write("%s%s\n" % ("END: Length of genomeList is ",len(genomeList)))

    ##### Other processing #####

    elif match_geneticCode:
        value = match_geneticCode.group(1)
        if value != '':
            geneticCode = value

    elif match_translateOnly:
        value = match_translateOnly.group(1)
        if value.lower() == 'yes' or value.lower() == 'true' or value.lower() == 'on':
            TRANSLATE_ONLY = True
        elif value.lower() == 'no' or value.lower() == 'false' or value.lower() == 'off' or value == '':
            TRANSLATE_ONLY = False
        else:
            if PHATE_WARNINGS == 'True':
                print("WARNING:  Invalid string following translate_only parameter in config file:", value)
            LOGFILE.write("%s%s\n" % ("Invalid string following translate_only parameter in config file: ", value))

    elif match_contig:
        value = match_contig.group(1)
        contigName = value


    ##### Gene Calls #####

    elif match_geneCaller:
        value = match_geneCaller.group(1)
        if value.lower() == 'phanotate':
            geneCaller = 'phanotate'
            CONSENSUS_CALLS_FILE = 'phanotate.cgc'
        elif value.lower() == 'consensus':
            geneCaller = 'consensus'
            CONSENSUS_CALLS_FILE = 'consensus.cgc'
        elif value.lower() == 'genemarks' or value.lower() == 'genemark':
            geneCaller = 'genemarks'
            CONSENSUS_CALLS_FILE = 'genemark.cgc'
        elif value.lower() == 'glimmer2':
            geneCaller = 'glimmer2'
            CONSENSUS_CALLS_FILE = 'glimmer.cgc'
        elif value.lower() == 'glimmer3' or value.lower() == 'glimmer':
            geneCaller = 'glimmer3'
            CONSENSUS_CALLS_FILE = 'glimmer.cgc'
        elif value.lower() == 'prodigal':
            geneCaller = 'prodigal'
            CONSENSUS_CALLS_FILE = 'prodigal.cgc'
        elif value.lower() == 'rast':
            geneCaller = 'rast'
            CONSENSUS_CALLS_FILE = 'rast.cgc'

    elif match_genemarksCalls:
        value = match_genemarksCalls.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
            genemarksCalls = True
        else:
            genemarksCalls = False

    elif match_prodigalCalls:
        value = match_prodigalCalls.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
            prodigalCalls = True
        else:
            prodigalCalls = False

    elif match_glimmerCalls:
        value = match_glimmerCalls.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
            glimmerCalls = True
        else:
            glimmerCalls = False

    elif match_phanotateCalls:
        value = match_phanotateCalls.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
            phanotateCalls = True
        else:
            phanotateCalls = False

    ##### BLAST #####

    elif match_blastpIdentity:
        value = match_blastpIdentity.group(1)
        if int(value) > int(MIN_BLASTP_IDENTITY) and int(value) <= 100:
            blastpIdentity = value

    elif match_blastpHitCount:
        value = match_blastpHitCount.group(1)
        if int(value) > 0 and int(value) <= int(MAX_BLASTP_HIT_COUNT):
            blastpHitCount = value

    elif match_blastnHitCount:
        value = match_blastnHitCount.group(1)
        if int(value) > 0 and int(value) <= int(MAX_BLASTN_HIT_COUNT):
            blastnHitCount = value

    elif match_ncbiVirusBlast:
        value = match_ncbiVirusBlast.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
            ncbiVirusBlast = True
        else:
            ncbiVirusBlast = False

    elif match_ncbiVirusProteinBlast:
        value = match_ncbiVirusProteinBlast.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
            ncbiVirusProteinBlast = True
        else:
            ncbiVirusProteinBlast = False

    elif match_keggVirusBlast:
        value = match_keggVirusBlast.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
             keggVirusBlast = True
        else:
             keggVirusBlast = False

    elif match_nrBlast:
        value = match_nrBlast.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
             nrBlast = True
        else:
             nrBlast = False 

    elif match_refseqProteinBlast:
        value = match_refseqProteinBlast.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
            refseqProteinBlast = True
        else:
            refseqProteinBlast = False

    elif match_refseqGeneBlast:
        value = match_refseqGeneBlast.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
            refseqGeneBlast = True
        else:
            refseqGeneBlast = False

    elif match_phantomeBlast:
        value = match_phantomeBlast.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
            phantomeBlast = True
        else:
            phantomeBlast = False 

    elif match_pvogsBlast:
        value = match_pvogsBlast.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
            pvogsBlast = True
        else:
            pvogsBlast = False

    elif match_uniparcBlast:
        value = match_uniparcBlast.group(1).lower()
        if value == 'true' or value == 'yes' or value == 'on':
            uniparcBlast = True
        else:
            uniparcBlast = False

    elif match_uniprotBlast:
        value = match_uniprotBlast.group(1).lower()
        if value == 'true' or value == 'yes' or value == 'on':
            uniprotBlast = True
        else:
            uniprotBlast = False

    elif match_swissprotBlast:
        value = match_swissprotBlast.group(1).lower()
        if value == 'true' or value == 'yes' or value == 'on':
            swissprotBlast = True
        else:
            swissprotBlast = False

    ##### HMM #####

    elif match_hmmProgram:
        value = match_hmmProgram.group(1)
        if value.lower() == 'jackhmmer':
            hmmProgram = 'jackhmmer'
        else:
            if PHATE_WARNINGS == 'True':
                print("WARNING: currenly only jackhmmer hmm search is supported; running jackhmmer")
            hmmProgram = HMM_PROGRAM_DEFAULT 

    elif match_ncbiVirusHmm:
        value = match_ncbiVirusHmm.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
            ncbiVirusHmm = True
        else:
            ncbiVirusHmm = False

    elif match_ncbiVirusProteinHmm:
        value = match_ncbiVirusProteinHmm.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
            ncbiVirusProteinHmm = True
        else:
            ncbiVirusProteinHmm = False

    elif match_keggVirusHmm:
        value = match_keggVirusHmm.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
             keggVirusHmm = True
        else:
             keggVirusHmm = False

    elif match_nrHmm:
        value = match_nrHmm.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
             nrHmm = True
        else:
             nrHmm = False 

    elif match_refseqProteinHmm:
        value = match_refseqProteinHmm.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
            refseqProteinHmm = True
        else:
            refseqProteinHmm = False

    elif match_phantomeHmm:
        value = match_phantomeHmm.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
            phantomeHmm = True
        else:
            phantomeHmm = False 

    elif match_pvogsHmm:
        value = match_pvogsHmm.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
            pvogsHmm = True
        else:
            pvogsHmm = False

    elif match_uniparcHmm:
        value = match_uniparcHmm.group(1).lower()
        if value == 'true' or value == 'yes' or value == 'on':
            uniparcHmm = True
        else:
            uniparcHmm = False

    elif match_uniprotHmm:
        value = match_uniprotHmm.group(1).lower()
        if value == 'true' or value == 'yes' or value == 'on':
            uniprotHmm = True
        else:
            uniprotHmm = False

    elif match_swissprotHmm:
        value = match_swissprotHmm.group(1).lower()
        if value == 'true' or value == 'yes' or value == 'on':
            swissprotHmm = True
        else:
            swissprotHmm = False

    elif match_refseqGeneHmm:
        pass  # Not yet in service

    ##### PSAT #####

    elif match_psatAnnotation:
        value = match_psatAnnotation.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
            psatAnnotation = True
        else:
            psatAnnotation = False

    elif match_psatFile:
        value = match_psatFile.group(1)
        if value != '':
            PSAT_FILE = value 
            PSAT = True   # Yes, a psat file will be passed to subordinate code
        else:
            PSAT_FILE = ""
            if PHATE_WARNINGS:
                print("multiPhate says, WARNING:  PSAT_FILE is ",PSAT_FILE)
        #LOG.write("%s%s\n" % ("PSAT_FILE is ",PSAT_FILE))

    ##### DEPENDENT CODE LOCATIONS #####

    elif match_blastPlusHome:
        if match_blastPlusHome.group(1) != '':
            os.environ["BLAST_HOME"] = match_blastPlusHome.group(1)

    elif match_embossHome:
        if match_embossHome.group(1) != '':
            os.environ["EMBOSS_PHATE_HOME"] = match_embossHome.group(1) 

    elif match_tRNAscanSEhome:
        if match_tRNAscanSEhome.group(1) != '':
            os.environ["tRNAscanSE_HOME"] = match_tRNAscanSEhome.group(1)

    elif match_glimmerHome:
        if match_glimmerHome.group(1) != '':
            os.environ["GLIMMER_PATH"] = match_glimmerHome.group(1)
 
    elif match_prodigalHome:
        if match_prodigalHome.group(1) != '':
            os.environ["PRODIGAL_PATH"] = match_prodigalHome.group(1)

    elif match_phanotateHome:
        if match_phanotateHome.group(1) != '':
            os.environ["PHANOTATE_PATH"] = match_phanotateHome.group(1)

    elif match_genemarkHome:
        if match_genemarkHome.group(1) != '':
            os.environ["GENEMARKS_PATH"] = match_genemarkHome.group(1)

    ##### DATABASE LOCATIONS #####

    elif match_ncbiVirusDatabase:
        if match_ncbiVirusDatabase.group(1) != '':
            os.environ["NCBI_VIRUS_BLAST_HOME"] = match_ncbiVirusDatabase.group(1) 

    elif match_refseqGeneDatabase:
        if match_refseqGeneDatabase.group(1) != '':
            os.environ["REFSEQ_GENE_BLAST_HOME"] = match_refseqGeneDatabase.group(1) 
  
    elif match_ncbiVirusProteinDatabase:
        if match_ncbiVirusProteinDatabase.group(1) != '':
            os.environ["NCBI_VIRUS_PROTEIN_BLAST_HOME"] = match_ncbiVirusProteinDatabase.group(1) 
  
    elif match_keggVirusDatabase:
        if match_keggVirusDatabase.group(1) != '':
            os.environ["KEGG_VIRUS_BLAST_HOME"] = match_keggVirusDatabase.group(1) 
            os.environ["KEGG_VIRUS_BASE_DIR"] = os.path.dirname(match_keggVirusDatabase.group(1)) + '/' 

    elif match_phantomeDatabase:
        if match_phantomeDatabase.group(1) != '':
            os.environ["PHANTOME_BLAST_HOME"] = match_phantomeDatabase.group(1) 
            os.environ["PHANTOME_BASE_DIR"] = os.path.dirname(match_phantomeDatabase.group(1)) + '/' 
 
    elif match_pvogsDatabase:
        if match_pvogsDatabase.group(1) != '':
            os.environ["PVOGS_BLAST_HOME"] = match_pvogsDatabase.group(1) 
            os.environ["PVOGS_HMM_HOME"]   = match_pvogsDatabase.group(1) 

    elif match_swissprotDatabase:
        if match_swissprotDatabase.group(1) != '':
            os.environ["SWISSPROT_BLAST_HOME"] = match_swissprotDatabase.group(1) 

    elif match_refseqProteinDatabase:
        if match_refseqProteinDatabase.group(1) != '':
            os.environ["REFSEQ_PROTEIN_BLAST_HOME"] = match_refseqProteinDatabase.group(1)

    elif match_nrDatabase:
        if match_nrDatabase.group(1) != '':
            os.environ["NR_BLAST_HOME"] = match_nrDatabase.group(1) 

    ##### VERBOSITY #####

    elif match_phateWarnings:
        value = match_phateWarnings.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
            os.environ["PHATE_WARNINGS"] = 'True' 
        else:
            os.environ["PHATE_WARNINGS"] = 'False'

    elif match_phateMessages:
        value = match_phateMessages.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
            os.environ["PHATE_MESSAGES"] = 'True' 
        else:
            os.environ["PHATE_MESSAGES"] = 'False' 

    elif match_phateProgress:
        value = match_phateProgress.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
            os.environ["PhATE_PROGRESS"] = 'True' 
        else:
            os.environ["PhATE_PROGRESS"] = 'False' 

    elif match_cgcWarnings:
        value = match_cgcWarnings.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
            os.environ["CGC_WARNINGS"] = 'True' 
        else:
            os.environ["CGC_WARNINGS"] = 'False' 
 
    elif match_cgcMessages:
        value = match_cgcMessages.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
            os.environ["CGC_MESSAGES"] = 'True' 
        else:
            os.environ["CGC_MESSAGES"] = 'False' 

    elif match_cgcProgress:
        value = match_cgcProgress.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
            os.environ["CGC_PROGRESS"] = 'True' 
        else:
            os.environ["CGC_PROGRESS"] = 'False' 

    elif match_cleanRawData:
        value = match_cleanRawData.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
            os.environ["CLEAN_RAW_DATA"] = 'True' 
        else:
            os.environ["CLEAN_RAW_DATA"] = 'False' 

    else:
        LOG.write("%s%s\n" % ("ERROR: Unrecognized line in config file: ", cLine))
        print("ERROR: unrecognized line in config file:", cLine)


LOG.write("%s\n" % ("Input parameters and configurables:"))
LOG.write("%s%s\n" % ("   BASE_DIR is ",os.environ["BASE_DIR"]))
LOG.write("%s%s\n" % ("   PHATE_BASE_DIR is ",os.environ["PHATE_BASE_DIR"]))
LOG.write("%s%s\n" % ("   DATABASE_DIR is ",os.environ["DATABASE_DIR"]))
LOG.write("%s%s\n" % ("   SOFTWARE_DIR is ",os.environ["SOFTWARE_DIR"]))
LOG.write("%s%s\n" % ("   GENE_FILE: ", GENE_FILE))
LOG.write("%s%s\n" % ("   PROTEIN_FILE: ", PROTEIN_FILE))

LOG.write("%s%s\n" % ("   geneticCode: ", geneticCode))
LOG.write("%s%s\n" % ("   Status of boolean TRANSLATE_ONLY is ",TRANSLATE_ONLY))
LOG.write("%s%s\n" % ("   geneCaller is ",geneCaller))
LOG.write("%s%s\n" % ("   genemarksCalls is ",genemarksCalls))
LOG.write("%s%s\n" % ("   prodigalCalls is ",prodigalCalls))
LOG.write("%s%s\n" % ("   glimmerCalls is ",glimmerCalls))
LOG.write("%s%s\n" % ("   phanotateCalls is ",phanotateCalls))
LOG.write("%s%s\n" % ("   CONSENSUS_CALLS_FILE is ",CONSENSUS_CALLS_FILE))

LOG.write("%s%s\n" % ("   blastpIdentity is ",blastpIdentity))
LOG.write("%s%s\n" % ("   blastpHitCount is ",blastpHitCount))
LOG.write("%s%s\n" % ("   blastnHitCount is ",blastnHitCount))
LOG.write("%s%s\n" % ("   ncbiVirusBlast is ",ncbiVirusBlast))
LOG.write("%s%s\n" % ("   ncbiVirusProteinBlast is ",ncbiVirusProteinBlast))
LOG.write("%s%s\n" % ("   keggVirusBlast is ",keggVirusBlast))
LOG.write("%s%s\n" % ("   nrBlast is ",nrBlast))
LOG.write("%s%s\n" % ("   refseqProteinBlast is ",refseqProteinBlast))
LOG.write("%s%s\n" % ("   refseqGeneBlast is ",refseqGeneBlast))
LOG.write("%s%s\n" % ("   phantomeBlast is ",phantomeBlast))
LOG.write("%s%s\n" % ("   pvogsBlast is ",pvogsBlast))
LOG.write("%s%s\n" % ("   uniparcBlast is ",uniparcBlast))
LOG.write("%s%s\n" % ("   swissprotBlast is ",swissprotBlast))

LOG.write("%s%s\n" % ("   hmmProgram is ",hmmProgram))
LOG.write("%s%s\n" % ("   ncbiVirusHmm is ",ncbiVirusHmm))
LOG.write("%s%s\n" % ("   ncbiVirusProteinHmm is ",ncbiVirusProteinHmm))
LOG.write("%s%s\n" % ("   keggVirusHmm is ",keggVirusHmm))
LOG.write("%s%s\n" % ("   nrHmm is ",nrHmm))
LOG.write("%s%s\n" % ("   refseqProteinHmm is ",refseqProteinHmm))
LOG.write("%s%s\n" % ("   refseqGeneHmm is ",refseqGeneHmm))
LOG.write("%s%s\n" % ("   phantomeHmm is ",phantomeHmm))
LOG.write("%s%s\n" % ("   pvogsHmm is ",pvogsHmm))
LOG.write("%s%s\n" % ("   uniparcHmm is ",uniparcHmm))
LOG.write("%s%s\n" % ("   swissprotHmm is ",swissprotHmm))

LOG.write("%s%s\n" % ("   blast+ home is ",os.environ["BLAST_HOME"]))
LOG.write("%s%s\n" % ("   emboss home is ",os.environ["EMBOSS_PHATE_HOME"]))
LOG.write("%s%s\n" % ("   tRNAscanSE home is ",os.environ["tRNAscanSE_HOME"]))
LOG.write("%s%s\n" % ("   glimmer home is ",os.environ["GLIMMER_PATH"]))
LOG.write("%s%s\n" % ("   prodigal home is ",os.environ["PRODIGAL_PATH"]))
LOG.write("%s%s\n" % ("   phanotate home is ",os.environ["PHANOTATE_PATH"]))
LOG.write("%s%s\n" % ("   genemark home is ",os.environ["GENEMARKS_PATH"]))

LOG.write("%s%s\n" % ("   ncbi virus database is located in ",os.environ["NCBI_VIRUS_BLAST_HOME"]))
LOG.write("%s%s\n" % ("   ncbi virus protein databases is located in ",os.environ["NCBI_VIRUS_PROTEIN_BLAST_HOME"]))
LOG.write("%s%s\n" % ("   refseq gene database is located in ",os.environ["REFSEQ_GENE_BLAST_HOME"]))
LOG.write("%s%s\n" % ("   kegg virus database is located in ",os.environ["KEGG_VIRUS_BLAST_HOME"]))
LOG.write("%s%s\n" % ("   phantome database is located in ",os.environ["PHANTOME_BLAST_HOME"]))
LOG.write("%s%s\n" % ("   pVOGs database is located in ",os.environ["PVOGS_BLAST_HOME"]))
LOG.write("%s%s\n" % ("   swissprot database is located in ",os.environ["SWISSPROT_BLAST_HOME"]))
LOG.write("%s%s\n" % ("   refseq protein database is located in ",os.environ["REFSEQ_PROTEIN_BLAST_HOME"]))
LOG.write("%s%s\n" % ("   NR database is located in ",os.environ["NR_BLAST_HOME"]))

LOG.write("%s%s\n" % ("   phate warnings is set to ",os.environ["PHATE_WARNINGS"]))
LOG.write("%s%s\n" % ("   phate messages is set to ",os.environ["PHATE_MESSAGES"]))
LOG.write("%s%s\n" % ("   phate progress is set to ",os.environ["PHATE_PROGRESS"]))
LOG.write("%s%s\n" % ("   cgc warnings is set to ",os.environ["CGC_WARNINGS"]))
LOG.write("%s%s\n" % ("   cgc messages is set to ",os.environ["CGC_MESSAGES"]))
LOG.write("%s%s\n" % ("   cgc progress is set to ",os.environ["CGC_PROGRESS"]))
LOG.write("%s%s\n" % ("   clean raw data is set to ",os.environ["CLEAN_RAW_DATA"]))

LOG.write("%s%s\n" % ("   PSAT is ",PSAT))
LOG.write("%s%s\n" % ("   PSAT_FILE is ",PSAT_FILE))

LOG.write("%s%s\n" % ("Number of genomes to be processed: ",len(genomeList)))
LOG.write("%s\n" % ("List of genomes to be processed:"))
for genome in genomeList:
    LOG.write("%s%c%s%c%s%c%s%c%s%c%s\n" % (genome["genomeNumber"],' ',genome["genomeName"],' ',genome["genomeType"],' ',genome["genomeSpecies"],' ',genome["genomeFile"],' ',genome["outputSubdir"]))

##### BEGIN MAIN ########################################################################################

# For each genome, create a phate.config file for running phate_runPipeline.py

nextConfigFile = ""
configList = []  # List of config filenames
for genome in genomeList:
    match_root = re.search(p_root,genome["genomeFile"])
    if match_root:
        genomeRoot = match_root.group(1)
        nextConfigFile = genomeRoot + '.config'
        configList.append(nextConfigFile)
        genomeFile = PIPELINE_INPUT_DIR + genome["genomeFile"]
    else:
        if PHATE_WARNINGS:
            print("multiPhate says, WARNING: Fasta filename not recognized for genome", genome["genomeName"])
            print("   Expected fasta filename extension: .fasta")
        LOG.write("%s%s%s\n" % ("WARNING: problem with fasta file: ",genome["genomeFile"]," This genome not processed!"))
        continue
    NEXT_CONFIG = open(nextConfigFile,"w")
    NEXT_CONFIG.write("%s%s%s%s\n" % ("# Config file for genome ",genome["genomeNumber"]," AKA ",genome["genomeName"])) 
    NEXT_CONFIG.write("\n%s\n" % ("# Processing information")) 
    NEXT_CONFIG.write("%s%c%s%c\n" % ("genome_file=","'",genome["genomeFile"],"'"))
    NEXT_CONFIG.write("%s%c%s%c\n" % ("genome_type=","'",genome["genomeType"],"'"))
    NEXT_CONFIG.write("%s%c%s%c\n" % ("species=","'",genome["genomeSpecies"],"'"))
    NEXT_CONFIG.write("%s%c%s%c\n" % ("name=","'",genome["genomeName"],"'"))
    NEXT_CONFIG.write("%s%c%s%c\n" % ("output_subdir=","'",genome["outputSubdir"],"'")) 
    NEXT_CONFIG.write("%s%c%s%c\n" % ("genetic_code=","'",geneticCode,"'")) 
    NEXT_CONFIG.write("%s%c%s%c\n" % ("translate_only=","'",TRANSLATE_ONLY,"'")) 
    NEXT_CONFIG.write("\n%s\n" % ("# Gene callers")) 
    NEXT_CONFIG.write("%s%c%s%c\n" % ("gene_caller=","'",geneCaller,"'")) 
    NEXT_CONFIG.write("%s%c%s%c\n" % ("genemarks_calls=","'",genemarksCalls,"'")) 
    NEXT_CONFIG.write("%s%c%s%c\n" % ("prodigal_calls=","'",prodigalCalls,"'")) 
    NEXT_CONFIG.write("%s%c%s%c\n" % ("glimmer_calls=","'",glimmerCalls,"'")) 
    NEXT_CONFIG.write("%s%c%s%c\n" % ("phanotate_calls=","'",phanotateCalls,"'")) 
    NEXT_CONFIG.write("\n%s\n" % ("# Blast parameters")) 
    NEXT_CONFIG.write("%s%c%s%c\n" % ("blast_identity=","'",blastpIdentity,"'")) 
    NEXT_CONFIG.write("%s%c%s%c\n" % ("blastn_hit_count=","'",blastnHitCount,"'")) 
    NEXT_CONFIG.write("%s%c%s%c\n" % ("blastp_hit_count=","'",blastpHitCount,"'")) 
    NEXT_CONFIG.write("%s\n" % ("# Blast databases to set")) 
    NEXT_CONFIG.write("%s\n" % ("# Genomes: ncbi_virus_blast")) 
    NEXT_CONFIG.write("%s%c%s%c\n" % ("ncbi_virus_blast=","'",ncbiVirusBlast,"'")) 
    NEXT_CONFIG.write("%s\n" % ("# Genes: refseq_gene_blast")) 
    NEXT_CONFIG.write("%s%c%s%c\n" % ("refseq_gene_blast=","'",refseqGeneBlast,"'")) 
    NEXT_CONFIG.write("%s\n" % ("# Proteins: ncbi_virus_protein_blast, kegg_virus_blast, nr_blast, refseq_protein_blast, phantome_blast, pvogs_blast")) 
    NEXT_CONFIG.write("%s%c%s%c\n" % ("ncbi_virus_protein_blast=","'",ncbiVirusProteinBlast,"'")) 
    NEXT_CONFIG.write("%s%c%s%c\n" % ("kegg_virus_blast=","'",keggVirusBlast,"'")) 
    NEXT_CONFIG.write("%s%c%s%c\n" % ("phantome_blast=","'",phantomeBlast,"'")) 
    NEXT_CONFIG.write("%s%c%s%c\n" % ("pvogs_blast=","'",pvogsBlast,"'")) 
    NEXT_CONFIG.write("%s%c%s%c\n" % ("swissprot_blast=","'",swissprotBlast,"'")) 
    NEXT_CONFIG.write("%s%c%s%c\n" % ("refseq_protein_blast=","'",refseqProteinBlast,"'")) 
    NEXT_CONFIG.write("%s%c%s%c\n" % ("nr_blast=","'",nrBlast,"'")) 
    NEXT_CONFIG.write("\n%s\n" % ("# HMM database(s) to set")) 
    NEXT_CONFIG.write("%s%c%s%c\n" % ("hmm_program=","'",hmmProgram,"'")) 
    NEXT_CONFIG.write("%s%c%s%c\n" % ("pvogs_hmm=","'",pvogsHmm,"'")) 
    NEXT_CONFIG.write("\n%s\n" % ("# Other, external annotation")) 
    NEXT_CONFIG.write("%s\n" % ("# If you have PSAT output, include the file")) 
    NEXT_CONFIG.write("%s\n" % ("# Otherwise, set to 'false' and set psat file to ''")) 
    NEXT_CONFIG.write("%s%c%s%c\n" % ("psat_annotation=","'",psatAnnotation,"'")) 
    NEXT_CONFIG.write("%s%c%s%c\n" % ("psat_file=","'",PSAT_FILE,"'")) 
    NEXT_CONFIG.close() 

# Run the pipeline (phate_runPipeline.py) over each genome; The code below runs in serial; Modify this section to implement in parallel on your cluster

LOG.write("%s%s\n" % ("Processing genomes through PhATE. Begin processing at ",datetime.datetime.now()))
for configFile in configList:
    LOG.write("%s%s\n" % ("Running PhATE using genome config file ",configFile))
    command = "python " + PHATE_BASE_DIR + PHATE_PIPELINE_CODE + " " + configFile 
    LOG.write("%s%s\n" % ("Command is ",command))
    LOG.write("%s%s\n" % ("Begin processing at ",datetime.datetime.now()))
    result = os.system(command)
    LOG.write("%s%s\n" % ("End processing at ",datetime.datetime.now()))

##### CLEAN UP

LOG.write("%s%s\n" % ("End log file ",datetime.datetime.now()))
LOG.close()
