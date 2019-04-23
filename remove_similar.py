import argparse
import itertools
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


#input_file - name of input file in fasta format
#Removes sequences with p-distance less than cutoff and more than cutoff1
def remove_sim_seq(input_file, cutoff, cutoff1):

    #cutoff = float(input('Input % distance cut-off > '))
    keys_alignment, alignment = parse_input_file(input_file)
    print(keys_alignment[0])
    print("Total sequence number: %d" % len(keys_alignment))

    for i, key in enumerate(keys_alignment):
        print("Progress iteration %d, left elements: %d" % (i, len(alignment)))
        if key in alignment:
            for k in itertools.islice(keys_alignment, i + 1, None):
                if k in alignment:
                    comp = compare_seq(alignment[key], alignment[k])
                    if comp < cutoff or comp > cutoff1:
                        print(comp, key, k)
                        del alignment[k]

    new_fn = "%s_%s.fasta" % (input_file.replace(".fasta", ""), cutoff)
    with open(new_fn, 'w') as f:
        for key in keys_alignment:
            if key in alignment:
                f.write(key)
                f.write("\n")
                f.write(alignment[key])
                f.write("\n")
    f.close()            
    return(new_fn)

#s1,s2 - strings
#returns similarity between 2 sequences BUT these sequences can be not aligned
def compare_seq(s1, s2):
    #number of equal positions
    eq = 0
    #number of unequal positions
    neq = 0
    s1 = s1.lower()
    s2 = s2.lower()
    
    if len(s1) > len(s2):
        temp = s2
        s2 = s1
        s1 = temp

    for n, s1_char in enumerate(s1):
        #we do not compare columns with only gaps
        if s1_char == "-" or s2[n] == "-":
            continue
            
        if s1_char == s2[n]:
            eq += 1
        else:
            neq += 1
        #print(s1_char, s2[n])
        #print(neq, eq)

    if neq == 0 and eq == 0: #sequences consist of gaps
        return 100
    else:
        return 100 * float(neq) / (neq + eq)
#input_file - name of input fasta-file
#returns sequence names as a list 'keys' and dictionary 'alignment': alignment[key] = ''
def parse_input_file(input_file):
    
    alignment = {}
    with open(input_file, 'r') as f:
        keys = []
        for line in f:
            temp = ''
            k=0
            if line.startswith('>'):
                key = line.strip('\n')
                if key in keys:
                    k=1
                    continue
                keys.append(key)
                alignment[key] = ""
            else: 
                if k==1:
                    k=0
                    continue
                else:
                    alignment[key] += line.splitlines()[0]
    #print(alignment[keys[0]])
    return keys, alignment
        
        
        



#finds ORF in sequence
#input - nucleotide sequence (type string)
#output - nucleotide sequence (only found ORF) (type string)
def find_ORF(seq):
    #threshold - minimal length of putative ORF
    if int(len(seq)/4)> 300:
        threshold = int(len(seq)/3.5)
    else: threshold = 300


    
    stop_codons = ['TAG', 'TAA', 'TGA', 'tag', 'taa', 'tga']
    #template seq
    seq1=seq[:]
    
    for nc in seq1:
        #counter of nucleotides
        i=0
        #equal to 1 if stop codon found
        k=0
        while i+3 <len(seq):
            #is it a stop codon?
            if str(seq[i:i+3]) in stop_codons:
                k = 1
                
                #if the length of seq is bigger than threshold ORF is found
                if i > threshold:
                    if len(seq) > 270:
                        return seq
                    else: 
                        return None
                #if the length of seq before stop codon is too short
                #the open reading frame might be shifted
                #we delete the first nucleotide to artificially shift it
                else:
                    seq=seq[1:]
                    break
            i+=3
        if k==0:
            if len(seq) > 270:
                return seq
            else:
                return None


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-input", "--input_file", type=str,
                        help="Input file in fasta format")
    parser.add_argument("-min", "--min_distance", type=int,
                        help="Minimal pairwise distance between sequences.\
                        If p-distance is lower than min_distance sequence \
                        with higher serial number will be removed from the dataset")
    parser.add_argument("-max", "--max_distance", type=int,
                        help="Maximal pairwise distance between sequences.\
                        If p-distance is higher than max_distance sequence \
                        with higher serial number will be removed from the dataset")
                            
    args = parser.parse_args()

    if not len(sys.argv) == 4:
        print('Please use "python remove_similar --help" for help ')

    else:
        remove_sim_seq(args.input_file, args.min_distance, args.max_distance)