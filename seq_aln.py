#!/usr/bin/env python

from Bio import SeqIO
from Bio.Align.Applications import ClustalwCommandline
from Bio.Seq import Seq
import sys, os
from Bio.SeqRecord import SeqRecord
import subprocess

def remove_stop_codons(dna_seq):#use with caution!

    '''removes typical stop codons from the ends of dna sequences. Takes a
    dna sequence as a string and returns a modified dna seq as a string'''
    stop_codons=["TAA","TAG","TGA"] #standard genetic code
    if dna_seq[-3:] in stop_codons:
        dna_seq=dna_seq[:-3]
    return dna_seq

def file_list(pattern=""):
    '''returns a list of filenames based on a pattern'''
    #os.system("ls " + pattern + "> tmp")
    cmd = "ls " + pattern + "> tmp; exit 0"
    subprocess.check_output(cmd,shell=True,stderr=subprocess.STDOUT)
    file_list = [name.strip() for name in open("tmp")]
    return file_list

#this is a comment
def split_dna(dna_seq):

    '''breaks a dna seq string into a list of triplets (codons)
    ["atg","taa"]'''
    codons = []
    triplet = ""
    for nucleotide in dna_seq:
        triplet =triplet + nucleotide
        if len(triplet)==3:
            codons.append(triplet)
            triplet = ""
    return codons

def check_seq(pep_dic,dna_dic):
    '''compares a conceptual translation to the peptide sequence in the aligment
    to make sure that the dna sequence given actually encodes the protein under
    the same name'''
    translation_dic = {}
    #make a dic of conceptually translated sequences
    missing_list = []
    #crashing here because of gaps?
    for dna_names in dna_dic: 
        translation_dic[dna_names]=str(dna_dic[dna_names].seq.translate())
    for pep_names in pep_dic:
        errors=[]
        error_file ="error_log.txt"
        pep_seq = pep_dic[pep_names].seq.tostring().replace("-","")
        #check for name miss-matches and report then crash if found
        try:
            dna_translation = translation_dic[pep_names]
        #catch sequence protein seq that dont have a DNA partner
        except KeyError:
            errors.append(pep_names)
        if len(errors)>0:
            error_handle = open(error_file,"w")
            errors_message = "".join(errors)
            error_handle.write(errors_message)
            error_handle.close()
            print "DNA and protein seq names don't match, check %s"%(error_file)
            print "The program will now 'crash'!"
            exit() #crash the program
        #check that the two sequences are the same!
        if pep_seq.upper().strip("*") != dna_translation.strip("*"):
            print pep_seq,dna_translation
            missing_list.append(pep_names)
    if len(missing_list) == 0:
        return None #success!
    else:
        return missing_list #problems with AA sequence itself
def find_gap(AA_seq):
    '''make a list of gap positions'''
    counter = 0
    position_list = []
    for aa in AA_seq:
        if aa == "-":
            position_list.append(counter)
        counter+=1
    return position_list

def insert_gaps(aa_seq,dna_seq):
    '''split dna into codons then insert gaps based on gap list. By using
    codons 1 AA postion = 1 DNA (codon position) and you can map the sequences
    to each other. Note that no codon/AA matching is done, that is why I pass the
    data to the check_seq function first'''
    codons=split_dna(dna_seq)
    gap_list = find_gap(aa_seq)
    for numbers in gap_list:
        codons.insert(numbers,"---")
    edited_dna = "".join(codons)
    gapped_dna = remove_stop_codons(edited_dna)
    return gapped_dna

def gap_a_lator(prot_aln_dic, dna_dic, outfile_name):

    '''takes protein alignments as a dictionary and their corresponding
    open reading frame DNA (with the same name!) and aligns the DNA based
    on the protein sequence alignment. This will perform checks and report
    errors if anything funny is going on.
    '''
    missing_list = check_seq(prot_aln_dic,dna_dic)
    if missing_list != None:
        print "problems found in:",missing_list
    for name in prot_aln_dic:
        pep_seq = prot_aln_dic[name].seq.tostring().upper()
        dna_seq = dna_dic[name].seq.tostring().upper()
        gapped_dna = insert_gaps(pep_seq,dna_seq)
        dna_dic[name].seq = Seq(gapped_dna)
    #outfile_name = "gapped.fna"
    outfile_handle = open(outfile_name,"w")
    for names in dna_dic:
        outfile_handle.write(dna_dic[names].format("fasta"))
    outfile_handle.close()

def make_protein_record(nuc_record):
    """Returns a new SeqRecord with the translated sequence (default table).
    proteins = (make_protein_record(nuc_rec) for nuc_rec in \
            SeqIO.parse("coding_sequences.fasta", "fasta"))
    SeqIO.write(proteins, "translations.fasta", "fasta")"""

    return SeqRecord(seq = nuc_record.seq.translate(cds=True), \
                     id = "" + nuc_record.id, \
                     description = "translation of CDS, using default table")

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
#collect gene names for all succesful alignments
good_genes = []

log = open('log.txt', 'w')
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

    info= "%s\t%s targets from %s strains\n"%(gene, len(target_seqs),
            len(strains))
    log.write(info)
    #this didnt work as never enough genes
    #if len(target_seqs) == len(strains)+1:
    good_genes.append(gene)
    tmp = open('./unaln/%s.faa'%gene, 'w')
    for rec in target_seqs:
        SeqIO.write(rec, tmp, 'fasta')
    tmp.close()
log.close()
#get read for algn
aln_files = file_list('./unaln/*.faa')
clustal_bin = '/usr/bin/clustalw'

for fil in aln_files:
    infile = fil
    outfile = "./aln/%s.aln.faa"%(fil.split('/')[-1].split('.')[0])#just name

    cline = ClustalwCommandline(clustal_bin, infile=infile,
            outfile=outfile, type="PROTEIN", output="FASTA",
            outorder="INPUT")

    #print cline
    #assert os.path.isfile(clustal_bin), "binary alignment path wrong"
    #stdout, stderr = cline()

#now get the DNA sequences
strain_dict_dna = {}
print "here1"
for fil in dna_files:
    spec = fil.split('./data/dna/')[1].split('_')[0]
    found = []
    #some genes are found twice only use first one
    tmp = open('tmp.fas', 'w')
    with open(fil) as f:
        for rec in SeqIO.parse(f, 'fasta'):
            if rec.id not in found:
                SeqIO.write(rec, tmp, 'fasta')
                found.append(rec.id)
    tmp.close()
    with open('tmp.fas') as f:
        tmp_dict = SeqIO.to_dict(SeqIO.parse(f, 'fasta'))
    strain_dict_dna[spec] = tmp_dict

#collect the DNA sequences for each gene in good_genes
strains_dna = strain_dict_dna.keys()

#get rid of reference
strains_dna.remove(reference)
ref_genes_dna = strain_dict_dna[reference]
good_genes_dna = []
target_seqs_dna = []

#I think this will break with missing species
print "here2"
for gene in good_genes:
    target_seqs_dna = [ref_genes_dna[gene]] #add the reference seq first
    for strain in strains_dna:
        key = "%s_%s"%(gene, strain)
        try:
            target = strain_dict_dna[strain][key]
            target_seqs_dna.append(target)
        except KeyError:
            pass
    #only do this if all strains are repressented
    #if len(target_seqs_dna) == len(strains_dna)+1:
    #should be == to good_genes
    good_genes_dna.append(gene)
    tmp = open('./unaln/%s.fas'%gene, 'w')
    for rec in target_seqs_dna:
        SeqIO.write(rec, tmp, 'fasta')
    tmp.close()

#check as pep from DNA should be equal length
#will fail because no genes conserved over all species
#assert good_genes == good_genes_dna

print "here3"
#now constrain dna alignment based on pep aln
#need seq dict for aln file
aln_files = file_list('./aln/*.aln.faa')
for aln in aln_files:
    dna_target = aln.replace('./aln','./unaln').replace('aln.faa','fas')
    pep_dict = {}
    for rec in SeqIO.parse(aln, 'fasta'):
        pep_dict[rec.id] = rec
    dna_dict = {}
    for rec in SeqIO.parse(dna_target, 'fasta'):
        dna_dict[rec.id] = rec
    
    outfile = './gapped_aln/'+ rec.id.split('_')[0] + "_gapped.aln.fas"
    #now do the alignment
    gap_a_lator(pep_dict, dna_dict, outfile)
