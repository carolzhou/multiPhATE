#!/usr/bin/env python

#############################################################################
#
# program: getDBs.py
#
# programmer:  C. E. Zhou
#
# Summary:  This script facilitates the downloading of databases to be used with multiPhATE.
#
# Most recent update:  18 March 2019
#
##############################################################################

import os, sys, re, time, datetime

##############################################################################
# CONSTANTS, BOOLEANS

BLAST               = False
NCBI_VIRUS_GENOME   = False
NCBI_REFSEQ_PROTEIN = False
NCBI_REFSEQ_GENE    = False
NCBI_SWISSPROT      = False
NR                  = False

# VARIABLES
blast               = ''
ncbi_virus_genome   = ''
ncbi_refseq_protein = ''
ncbi_refseq_gene    = ''
ncbi_swissprot      = ''
nr                  = ''
decision            = ''
blastPath           = '' # path to user's blast installation
cwd                 = '' # current working direcgtory

# Set up database directories
cwd              = os.getcwd()
dbDir            = cwd       + "/Databases/"
ncbiDir          = dbDir     + "NCBI/"
ncbiGenomeDir    = ncbiDir   + "Virus_Genome/"
ncbiProteinDir   = ncbiDir   + "Virus_Protein/"
nrDir            = dbDir     + "NR/"
refseqDir        = dbDir     + "Refseq/"
refseqProteinDir = refseqDir + "Protein/"
refseqGeneDir    = refseqDir + "Gene/"
swissprotDir     = dbDir     + "Swissprot/"

if not os.path.exists(dbDir):
    os.mkdir(dbDir)
if not os.path.exists(ncbiDir):
    os.mkdir(ncbiDir)
if not os.path.exists(ncbiGenomeDir):
    os.mkdir(ncbiGenomeDir)
if not os.path.exists(ncbiProteinDir):
    os.mkdir(ncbiProteinDir)
if not os.path.exists(refseqDir):
    os.mkdir(refseqDir)
if not os.path.exists(refseqProteinDir):
    os.mkdir(refseqProteinDir)
if not os.path.exists(refseqGeneDir):
    os.mkdir(refseqGeneDir)
if not os.path.exists(swissprotDir):
    os.mkdir(swissprotDir)

##############################################################################
# First, determine if user needs to download BLAST+.

print ("Let's download the databases you will need to run multiPhATE")
print ("First, you need blast+ in order to install several of the databases")
print ("Please confirm that you have downloaded and installed blast+: type 'y' (yes) or 'n' (no)")
blast = input()
if re.search('Y|y|Yes|yes|YES', blast):
    BLAST = True 
    print ("That's great! Now tell me where your blast executables are located.") 
    print ("If you installed them within your Conda environment, then you should find")
    print ("them in a path something like, \"/Users/myName/miniconda3/envs/multiPhate/bin/\"")
    print ("That directory should contain executables: blastn, blastp, blastx, and makeblastdb")
    print ("Please input the fully qualified path to blast+ : ")
    blastPath = input()
    # CHECK THAT THIS IS CORRECT !!!"
elif re.search('N|n|No|no|NO', blast):
    BLAST = False 
    print ("Please consult the README file for how to acquire and install BLAST+.") 
    print ("Note: The easiest way to install BLAST+ is within a Conda environment.")
    print ("Without blast+, we can still download some of the databases.")
    print ("Shall we continue? please respond 'y' or 'n': ")
    toContinue = input()
    if re.search('Y','y','yes','Yes','YES', toContinue):
        pass 
    else:
        print ("Bye")
        exit()
else:
    print ("That was not a correct response; please try again")
    exit()

##############################################################################
# Next, determine which databases to download

print ("Please indicate which databases you would like to download...")
print ("For each, type 'y' or 'n'")

#####
print ("NCBI Virus Genome database: ('y'/'n')")
ncbi_virus_genome = input()
if re.search('Y|y|yes|Yes|YES',ncbi_virus_genome):
    print ("Great, let's download NCBI Virus Genome")
    NCBI_VIRUS_GENOME = True
elif re.search('N|n|no|No|NO',ncbi_virus_genome):
    print ("Ok, we'll skip that one")
else:
    print ("That was not a correct response; please try again: type 'y' or 'n'")

#####
print ("NCBI Refseq Protein database: ('y'/'n')")
ncbi_refseq_protein = input()
if re.search('Y|y|yes|Yes|YES',ncbi_refseq_protein):
    print ("Great, let's download NCBI Virus Genome")
    NCBI_REFSEQ_PROTEIN = True
elif re.search('N|n|no|No|NO',ncbi_refseq_protein):
    print ("Ok, we'll skip that one")
else:
    print ("That was not a correct response; please try again: type 'y' or 'n'")

#####
print ("NCBI Refseq Gene database: ('y'/'n')")
ncbi_refseq_gene = input()
if re.search('Y|y|yes|Yes|YES',ncbi_refseq_gene):
    print ("Great, let's download NCBI Virus Genome")
    NCBI_REFSEQ_GENE = True
elif re.search('N|n|no|No|NO',ncbi_refseq_gene):
    print ("Ok, we'll skip that one")
else:
    print ("That was not a correct response; please try again: type 'y' or 'n'")

#####
print ("NCBI Swissprot database: ('y'/'n')")
ncbi_swissprot = input()
if re.search('Y|y|yes|Yes|YES',ncbi_swissprot):
    print ("Great, let's download NCBI Virus Genome")
    NCBI_SWISSPROT = True
elif re.search('N|n|no|No|NO',ncbi_swissprot):
    print ("Ok, we'll skip that one")
else:
    print ("That was not a correct response; please try again: type 'y' or 'n'")

#####
print ("Caution: the NR database is extremely large. ")
print ("If you already have it on disk, you are advised not to download here. ")
print ("Downloading NR will take a long time. ")
print ("NR database: ('y'/'n')")
nr = input()
if re.search('Y|y|yes|Yes|YES',nr):
    print ("Great, let's download NCBI Virus Genome")
    NR = True
elif re.search('N|n|no|No|NO',nr):
    print ("Ok, we'll skip that one")
else:
    print ("That was not a correct response; please try again: type 'y' or 'n'")

#####
if (BLAST or NCBI_VIRUS_GENOME or NCBI_REFSEQ_PROTEIN or NCBI_REFSEQ_GENE or NCBI_SWISSPROT or NR):
    print ("Ok, this is what we are going to download. ")
    if BLAST:
        print ("BLAST plus")
    if NCBI_VIRUS_GENOME:
        print ("NCBI Virus Genome database")
    if NCBI_REFSEQ_PROTEIN:
        print ("NCBI Refseq Protein database")
    if NCBI_REFSEQ_GENE:
        print ("NCBI Refseq Gene database")
    if NCBI_SWISSPROT:
        print ("NCBI Swissprot database")
    if NR:
        print ("NR database")
    print ("Print 'go' to proceed, or 'stop' to reconsider. ")
    print ("(You can always run this script again.) ")
    decision = input()
    if (re.search('go|Go|GO',decision)):
        print ("Ok, let's get started with downloading.")
        print ("Databases will be installed into the Databases/ folder within multiPhATE")
    else:
        print ("Ok, maybe some other time. Bye!")
        exit()
else:
    print ("You have selected no downloads. Have a happy day! :-)")
    exit()

##############################################################################
# Install BLAST+

# Install blast+ for user; NOT YET IN SERVICE
#if not BLAST:
#    command1 = "wget \"ftp//ftp.ncbi.nih.gov/blast/executables/blast+/${BLAST_VERSION}/ncbi- blast-${BLAST_VERSION}+-x64-linux.tar.gzi\""
#    command2 = "tar -zxf ncbi-blast-${BLAST_VERSION}+-x64-linux.tar.gz"
#    command3 = "cd ncbi-blast-${BLAST_VERSION}+/bin"
#    command4 = "pwd"

##############################################################################
# Install NCBI_VIRUS_GENOME

if NCBI_VIRUS_GENOME:
    os.chdir(ncbiGenomeDir)
    # VIRUS GENOME DB NOT FOUND - THIS ONE IS MANUAL 
    try:
        print ("Downloading NCBI Virus Genome database.")
        print ("This may take a while...")
        command = blastPath + "update_blastdb.pl" + ' ' + "virus_genome"
        #success = os.system(command)
        print ("NCBI Virus Genome database download complete.")
        print ("Formatting database for blast.")
        command = "makeblastdb -dbtype nucl -in viral.genomic"
        #success = os.system(command)
        print ("Database is formatted.")
    except BlastError:  
        print ("Command " + command + " failed; please check the location of your blast executables")
    os.chdir(cwd)

##############################################################################
# Install NCBI_REFSEQ_PROTEIN

if NCBI_REFSEQ_PROTEIN:
    os.chdir(ncbiProteinDir)
    try:
        print ("Downloading NCBI Refseq Protein database.")
        print ("This may take a while...")
        command = blastPath + "update_blastdb.pl" + ' ' + "refseq_protein"
        success = os.system(command)
        print ("NCBI Refseq Protein database download complete.")
        print ("Formatting database for blast.")
        command = "makeblastdb -dbtype prot -in refseq_protein"
        success = os.system(command)
        print ("Database is formatted.")
    except BlastError:  
        print ("Command " + command + " failed; please check the location of your blast executables")
    os.chdir(cwd)

##############################################################################
# Install NCBI_REFSEQ_GENE

if NCBI_REFSEQ_GENE:
    os.chdir(refseqGeneDir)
    try:
        print ("Downloading NCBI Refseq Gene database.")
        print ("This may take a while...")
        command = blastPath + "update_blastdb.pl" + ' ' + "refseqgene"
        success = os.system(command)
        print ("NCBI Refseq Gene database download complete.")
        print ("Formatting database for blast.")
        command = "makeblastdb -dbtype nucl -in refseqgene"
        success = os.system(command)
        print ("Database is formatted.")
    except BlastError:
        print ("Command " + command + " failed; please check the location of your blast executables")
    os.chdir(cwd)

##############################################################################
# Install SWISSPROT

if NCBI_SWISSPROT:
    os.chdir(swissprotDir)
    try:
        print ("Downloading Swissprot database.")
        print ("This may take a while...")
        command = blastPath + "update_blastdb.pl" + ' ' + "swissprot"
        success = os.system(command)
        print ("NCBI Swissprot database download complete.")
        print ("Formatting database for blast.")
        command = "makeblastdb -dbtype prot -in swissprot"
        success = os.system(command)
        print ("Database is formatted.")
    except BlastError:
        print ("Command " + command + " failed; please check the location of your blast executables")
    os.chdir(cwd)

##############################################################################
# Install NR

if NR:
    os.chdir(nrDir)
    try:
        print ("Downloading NCBI NR database.")
        print ("This may take a while...")
        command = blastPath + "update_blastdb.pl" + ' ' + "nr"
        success = os.system(command)
        print ("NCBI NR database download complete.")
        print ("Formatting database for blast.")
        command = "makeblastdb -dbtype nucl -in nr"
        success = os.system(command)
        print ("Database is formatted.")
    except BlastError:
        print ("Command " + command + " failed; please check the location of your blast executables")
    os.chdir(cwd)

##############################################################################

print ("Done!")

##############################################################################
##############################################################################

