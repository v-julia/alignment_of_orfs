import argparse
import os
from Bio import AlignIO


# Input - name of file with alignment in fasta-format
# Finds positions without gaps in reference sequence (the first in alignment)
# Removes the columns with gaps in reference sequence from alignment
# Saves new alignment to the file 'input_file_name_cut.fasta'

def cut_fasta(input_file):
    
    alignment = AlignIO.read(open(input_file), "fasta") # alignment object
    temp_seq = alignment[0].seq # template - reference sequence
    

    positions = [] # list of positions without gaps in reference sequence
    
    
    pos_st = 0 # start position of block without gaps in reference seq
    pos_end = 0 # end position of block without gaps in reference seq
    
    count = 0 #legnth of block 
    prev = '' #previous nucleotide

    #if k==1, we found block of gaps
    k = 0

    #searching blocks without gaps
    for nuc in temp_seq:
        if nuc=='-' and prev!='-':
            pos_end = count
            if pos_end !=0:
                positions.append([pos_st, pos_end])
                k = 1
        if nuc!='-' and prev=='-':
                pos_st = count
                k = 0
        count+=1
        prev = nuc
    if k == 0:
        positions.append([pos_st, len(temp_seq)])
    print(positions)
    
    #if no gaps in reference seq
    if len(positions) == 0:
        alignment1 = alignment
    #cutting regions without gaps
    else:
        #alignment1 = alignment[:,positions[0][0]:(positions[len(positions)-1][1])]
        alignment1 = alignment[:,positions[0][0]:(positions[0][1])]
        for i in range (1, len(positions)):
            alignment1 = alignment1 + alignment[:,positions[i][0]:(positions[i][1])]


    #If more than 10% of sequence are gaps, the sequence is deleted

    alignment_l = []
    for rec in alignment1:
        count_gap = rec.seq.count('-')
        #print(count_gap/len(rec.seq))
        if count_gap/len(rec.seq)<0.10:
            alignment_l.append(rec)
 
    alignment_new =  AlignIO.MultipleSeqAlignment(alignment_l)
    
    print('Number of sequences in alignment {}'.format(len(alignment_new)))
    #print( "Alignment length {0}".format(alignment_new.get_alignment_length()))

    out_file = input_file.replace('.fasta', '_cut.fasta')
    #out_file = '.'.join(input_file.split('.')[:-1]) + '_cut.fasta'
    AlignIO.write(alignment_new, open(out_file, 'w'), "fasta")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-input", "--input_file", type=str,
                        help="Input file")
    args = parser.parse_args()

    if not args.input_file:
        print("Please, use \"python cut_fasta.py --help\"")
    else:
        parse_gb(args.input_file)
