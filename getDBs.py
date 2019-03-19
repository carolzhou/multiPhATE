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

##############################################################################
# First, determine if user needs to download BLAST+.

print ("Let's download the databases you will need to run multiPhATE")
print ("First, tell me if you need to install blast+ (blast plus); type 'y' or 'n':")
print ("Download BLAST plus: ('y'/'n')")
blast = input()
if re.search('Y|y|Yes|yes|YES', blast):
    print ("Ok, then, let's install blast+.") 
    BLAST = True 
elif re.search('N|n|No|no|NO', blast):
    print ("Ok, then, we'll skip that one.") 
    BLAST = False 
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

else:
    print ("You have selected no downloads. Have a happy day! :-)")
    exit()

##############################################################################
# Install BLAST+

if BLAST:
    pass

##############################################################################
# Install NCBI_VIRUS_GENOME

if NCBI_VIRUS_GENOME:
    pass

##############################################################################
# Install NCBI_REFSEQ_PROTEIN

if NCBI_REFSEQ_PROTEIN:
    pass

##############################################################################
# Install NCBI_REFSEQ_GENE

if NCBI_REFSEQ_GENE:
    pass

##############################################################################
# Install SWISSPROT

if NCBI_SWISSPROT:
    pass

##############################################################################
# Install NR

if NR:
    pass

##############################################################################


##############################################################################
##############################################################################


