#!/usr/bin/env python

from Bio import SeqIO
from Bio.Align.Applications import ClustalwCommandline
import sys, os
from linux_control import file_list
from Bio.SeqRecord import SeqRecord

#this is a comment

#reference = 'S288C_reference_orf_trans_all.fasta'
def make_protein_record(nuc_record):
    """Returns a new SeqRecord with the translated sequence (default table).
    proteins = (make_protein_record(nuc_rec) for nuc_rec in \
            SeqIO.parse("coding_sequences.fasta", "fasta"))
    SeqIO.write(proteins, "translations.fasta", "fasta")"""

    return SeqRecord(seq = nuc_record.seq.translate(cds=True), \
                     id = "" + nuc_record.id, \
                     description = "translation of CDS, using default table")
'''
counter = 0
dna_files = file_list("./data/dna/*")
for fil in dna_files:
    outfile_n = './data/pep/%s'%(fil.
            split('./data/dna/')[1].split('.')[0]) + '.fsa'
    outfile = open(outfile_n, 'w')
    with open(fil) as f:
        found = []
        for rec in SeqIO.parse(f, 'fasta'):
            try:
                pep = make_protein_record(rec)
#sometimes they have the same name just grab one record
                if pep.id not in found:
                    SeqIO.write(pep, outfile, 'fasta')
                    found.append(pep.id)
            except:
                #print rec.id
                counter +=1
    print fil
    print counter

    #reset the counter
    counter = 0
    outfile.close()
'''
#make a dictionary to hold dictionaries of names
strain_dict = {}

pep_files = file_list("./data/pep/*")
for fil in pep_files:
    spec = fil.split('./data/pep/')[1].split('_')[0]
    with open(fil) as f:
        tmp_dict = SeqIO.to_dict(SeqIO.parse(f, 'fasta'))
    strain_dict[spec] = tmp_dict

strains = strain_dict.keys()
reference = 'S288C'
strains.remove(reference)

#get sequences ready for alignments
ref_genes = strain_dict[reference]
for gene in ref_genes:
    target_seqs = [ref_genes[gene]] #add the reference seq first
    for strain in strains:
        key = "%s_%s"%(gene, strain)
        try:
            target = strain_dict[strain][key]
            target_seqs.append(target)
        except KeyError:
            pass
    #only do this if all strains are repressented 
    if len(target_seqs) == len(strains)+1:
        tmp = open('./unaln/%s.faa'%gene, 'w')
        for rec in target_seqs:
            SeqIO.write(rec, tmp, 'fasta')
        tmp.close()

#get read for algn
aln_files = file_list('./unaln/*.faa')
clustal_bin = '/usr/bin/clustalw'
opts = " -type=PROTEIN -outputorder=INPUT"

for fil in aln_files:
    infile = fil
    outfile = "./aln/%s.aln.faa"%(fil.split('/')[-1].split('.')[0])#just name

    cline = ClustalwCommandline(clustal_bin, infile=infile,
            outfile=outfile, type="PROTEIN", output="FASTA",
            outorder="INPUT")

    #print cline
    #assert os.path.isfile(clustal_bin), "binary alignment path wrong"
    stdout, stderr = cline()


