#!/usr/bin/env python

from Bio import SeqIO
import sys
from linux_control import file_list

#this is a comment

#reference = 'S288C_reference_orf_trans_all.fasta'

target_files = file_list("./data/pep/*")

spec_list = [name.split('./data/pep/')[1].split('_')[0] for name in target_files]
'''
#check files available and make sure all species are there
target_species = ['S288C']
all_species = 1
for spec in target_species:
    try:
        assert spec in gene_list
    except AssertionError:
        print "missing species %s"%spec 
        all_species = None

if all_species == None:
    #get out of missing species TODO: fix this behaviour
    sys.exit()

#
'''
