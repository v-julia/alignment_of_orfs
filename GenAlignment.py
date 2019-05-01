import os
from os import system
from pathlib import Path
import sys
from Bio import SeqIO, Seq
from parser_gb import parse_gb
from remove_similar import remove_sim_seq
from remove_random import remove_random
from excise_target_region import find_target_region
from cut_fasta import cut_fasta

#This script:
#0) extracts parameters from 'config.txt' file stored in local directory
#1) converts genbank file to fasta format with sequences names in format "GenbankID_country_year"
#2) performs blast of sequences from genbank file against reference genome, saves sequences which overlap with 
#target genome (coordinates are defined in 'config.txt') are saved to file {input_fname}_exc.fasta
#3) adds reference sequences to {input_fname}_exc.fasta and saves the resulting dataset to {input_fname}_exc_wref.fasta
#4) multiple alignment using mafft > {input_fname}_exc_wref_aln.fasta
#5) Removes some sequences from the dataset in two ways:
#a) calculates p-distances between sequences pairs in loop, if p-distance < cutoff and p-distance > cutoff1
#removes the sequence with higher serial number. cutoff and cutoff1 must be written in config.txt
#Save sequences to {input_fname}_exc_wref_remove_cutoff.fasta
#b) divides all sequences into groups by the first 5 characters in GenBank ID, 
#randomly removes k% sequences in groups which size exceed m. k and m are defined by user.
#Saves sequences to {input_fname}_exc_wref_random_m_k.fasta
#6)cuts alignment by reference sequence, saves alignment to {input_fname}_exc_wref_random_m_k_cut.fasta

# Extracting parameters from config.txt

print('Reading config.txt ...')

config_file = open(os.path.join(sys.path[0],'config.txt'), 'r') #opens config.txt located in the same dir as script
for line in config_file:
    if line[0] == '#' or line == '\n':
        continue
    else:
        list_par = line.split('=')
        feature = list_par[0].strip(' ')
        value = list_par[1].strip(' ').strip('\n').strip('\'')
        if feature == 'input_file':
            input_file = '/'.join(value.split('\\')) #input file name
        if feature == 'min_length':
            min_length = int(value) #milinal length of sequence
        if feature == 'max_length':
            max_length = int(value) #maximal length of sequence
        if feature == 'reference':
            ref_fasta = '/'.join(value.split('\\')) #file with reference sequences
        if feature == "positions":
            rstart=int(value.split(",")[0]) #start of target region in reference sequence
            rend = int(value.split(",")[1]) #end of target region in reference sequence
        if feature == 'translation':
            translation = int(value) #amino-acid based alignment or not, in progress
        if feature == 'method':
            method = value #method of removing redundant sequences
            if (method != 'random' and method != 'similar'):
                print('Incorrect method')
        if feature == 'cutoff':
            cutoff = float(value) #minimal difference between 2 sequences if method=similar
        if feature == 'cutoff1':
            cutoff1 = float(value) #maximal difference between 2 sequences if method=similar
        if feature == 'path_to_blast': 
            path_to_blast = value #path to blast program
        if feature == 's':
            s = int(value) #intermediate files will be saved s==1
        if feature == 'path_to_mafft':
            path_to_mafft = value
config_file.close()
print('Done')


# 1) fasta from gb_file
print(input_file)
print('Parsing genbank file...')
# type - full capsid gene or short fragment
fasta0 = parse_gb(input_file, min_length, max_length) #name of file with annotated sequences extracted from gb file
print('Done')

#2) extracts sequences which contain target region from fasta0
fasta_targ = find_target_region(fasta0, ref_fasta, rstart, rend, path_to_blast)


#3) joins fasta-file with refseq_fasta
try:
    fasta_wref = fasta_targ.replace('.fasta', '_wref.fasta')

    ref = SeqIO.read(ref_fasta, "fasta")
    ref.seq = ref.seq[rstart-1:rend]
    ref_gb_id = ref.id.split('_')[0]

    records = list(SeqIO.parse(fasta_targ, "fasta"))

    for rec in records:             #removes duplicates of reference sequence
        if rec.id.split('_')[0] == ref_gb_id:
            records.remove(rec)
    SeqIO.write([ref]+records, fasta_wref, "fasta")
except:
    print('Error: Can\'t find file with  reference sequence ')

#4) aligns sequences

fasta_aln = fasta_wref.replace('.fasta', '_aln.fasta')

print("Aligning sequences")
if sys.platform == 'win32' or sys.platform == 'cygwin':
    system(('{} --op 60 --ep 9 --retree 1 ' + fasta_wref+ ' > ' + fasta_aln).format(path_to_mafft+'mafft.bat'))
else:
    system(('{} --op 60 --ep 9 --retree 1 ' + fasta_wref+ ' > ' + fasta_aln).format(path_to_mafft+'mafft'))

print("Done")


# 5a) Removes similar sequences from alignment to get rid of almost identical ones
if method == "similar":
    print('Removing similar seqs...')
    # cutoff - minimal p-distance between sequences
    # cutoff1 - maximal p-distance
    fasta_rem = remove_sim_seq(fasta_aln, cutoff, cutoff1)

# 5b) Divides sequences into groups by the first 6 characters in GB ID
#removes k% of sequences in groups with size > m
#k and m are defined by user
else:
    print("Removing random fraction of dataset")
    fasta_rem = remove_random(fasta_aln)
 
# 6)Cuts alignment based on reference sequence
#try:
print('Cutting sequences...')
cut_fasta(fasta_rem)
print('Done')
#except:
#    print('Error in cutting sequences')


# Deletes files generated on steps 1-5
if s==0:
    try:
        print('Deleting intermediate files...')

        out_dir = '/'.join(input_file.split('/')[:-1])+'/'

        Path(fasta0).unlink()
        Path(fasta_targ).unlink()
        Path(fasta_wref).unlink()
        Path(fasta_rem).unlink()
        Path(out_dir+'blast.out').unlink()
        Path(out_dir+'local_db.nhr').unlink()
        Path(out_dir+'local_db.nin').unlink()
        Path(out_dir+'local_db.nsq').unlink()

    except:
        print('Error: Can\'t  delete intermediate files')

print('The program\'s finished')
