import argparse
import pandas as pd
import os
import re
from Bio import GenBank
from matplotlib import pyplot as plt


def host_distr(input_file, output_dir, out_name):
    '''
    Plots distribution of values in host qualifier in genbank file
    '''

    host_dict = {}
    accession = re.compile(r"^ACCESSION\s+([a-z_A-Z0-9]+)") #accession number
    host = re.compile(r"^\s+\/host=\"([^\"]+)\"$") #qualifiers from FEATURE source field


    with open(input_file) as handle:
        for line in handle:

            m = host.match(line)
            if m:
                h = m.group(1)
                if h in host_dict.keys():
                    host_dict[h] += 1
                else:
                    host_dict[h] = 1
    print(host_dict)
    host_p = pd.DataFrame.from_dict(host_dict, orient='index',columns=['counts'])
    host_p = host_p.sort_index(ascending=True)
    host_p.to_csv(os.path.join(output_dir, args.output_name), header=True)
    host_p.plot(kind='bar')
        
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-input", "--input_file", type=str,
                        help="Input file in GenBank format",
                        required=True)
    parser.add_argument("-odir", "--output_dir", type=str,
                        help="Output directory to save the output file",
                        )
    parser.add_argument("-oname", "--output_name", type=str,
                        help="Name of output file",
                        required=True)
    args = parser.parse_args()
    if not args.output_dir:
        args.output_dir = ''
    host_distr(args.input_file, args.output_dir, args.output_name)
