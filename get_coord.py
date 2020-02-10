import argparse
import os
import re
from Bio import SeqIO
from Bio import GenBank



def get_utr_coord(input_file):
    '''
    Finds the coordinates of 5'NTR and 3'NTR in file with nucleotide sequences
    in GenBank format
    
    Input:
        input_file - path to file with nt seqs in GenBank format
    Outtput:
        writes coordinates to a new file in the same directory as input_file
    '''
    
    
    out_file_name = os.path.splitext(input_file)[0] + '_coord.txt'
    out_file = open(out_file_name, 'w')

    out_file.write('acc,st_5utr,e_5utr,st_cds,e_cds,st_3utr,e_3utr\n')

    REG_JOIN = re.compile(r'[\d]+\.\.[\d]+')
    COD_START = re.compile(r'codon_start=[123]')
    with open(input_file) as handle:
        for seq in GenBank.parse(handle):

            # string that will be written to file
            st = ''
            # accession number of seq
            acc = seq.accession[0]
            # length of sequence
            seq_length = int(seq.size)
            start_5utr = 0
            end_5utr = 0
            start_c = 0
            end_c = 0
            start_3utr = 0
            end_3utr = 0

            cds_features = list()
            for feature in seq.features:
                '''
                добавить 5'UTR, 3'UTR
                '''
                if feature.key == 'CDS':
                    cds_features.append(feature)
            if cds_features:
            
                cds_0_str = str(cds_features[0])
                cs = COD_START.search(cds_0_str)
                if cs:
                    codon_start = int(cs.group()[-1]) - 1
                else:
                    codon_start = 0
                
                
                if cds_features[0].location.startswith('join'):
                    loc = REG_JOIN.findall(cds_features[0].location)
                    st1, _ = loc[0].split('..')
                else:
                    st1, _ = cds_features[0].location.split('..')
                st1 = int(st1.strip('<').strip('>'))
                
                
                
                if cds_features[-1].location.startswith('join'):
                    loc = REG_JOIN.findall(cds_features[-1].location)
                    _, e2 = loc[1].split('..')
                else:
                    _, e2 = cds_features[-1].location.split('..')
                e2 = int(e2.strip('<').strip('>'))

                start_c = st1 + codon_start
                end_c = e2
                if st1 > 1:
                    start_5utr = 1
                    end_5utr = st1 - 1 + codon_start
                elif codon_start > 1:
                    start_5utr = 1
                    end_5utr = st1 - 1 + codon_start
                if e2 < seq_length:
                    start_3utr = e2 + 1
                    end_3utr = seq_length
                l_coord = [str(x) for x in [acc,
                                            start_5utr,
                                            end_5utr,
                                            start_c,
                                            end_c,
                                            start_3utr,
                                            end_3utr]]
                st = ','.join(l_coord) + '\n'
                out_file.write(st)
            else:
                print('No CDS for {}'.format(acc))

    out_file.close()
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-input", "--input_file", type=str,
                        help="Input file", required=True)
    args = parser.parse_args()

    get_utr_coord(args.input_file)