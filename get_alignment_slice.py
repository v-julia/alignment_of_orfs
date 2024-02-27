import argparse
import os
from Bio import AlignIO


def get_slice(input_file, start, end):
    '''
    Gets slice of alignment and saves it to fasta file
    
    input_file - name of file with alignment in fasta format
    start - start position of slice (starts from 1)
    end - end position of slice
    '''

    alignment = AlignIO.read(open(input_file), "fasta") # alignment object
    
    alignment_slice = alignment[:,start-1:end]
    out_file = str(os.path.splitext(input_file)[0] +'_' + str(start) + '-' + str(end) + '.fasta')
    AlignIO.write(alignment_slice, open(out_file, 'w'), "fasta")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-input", "--input_file", type=str,
                        help="Input file", required=True)
    parser.add_argument("-s", "--start", type=int,
                        help="Start of slice", required=True)
    parser.add_argument("-e", "--end", type=int,
                        help="End of slice", required=True)
    args = parser.parse_args()


    get_slice(args.input_file, args.start, args.end)
