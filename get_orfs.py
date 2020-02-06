import argparse
import copy
import os
import pandas as pd
import re
from Bio import SeqIO
from parser_gb import read_csv

def orf_coord(input_file, orf_map):
    '''
    Retrieves the coordinates of ORFs
    
    Input:
        input_file - file with nucleotide sequences in genbank-format
        orf_map - csv file annotation of orfs and their codes
    Output:
        coord_file - file with coordinates
    '''

    # dictionary with annotations of ORFs
    orf_dict = read_csv(orf_map)
    # possible ORFs
    orf_types = ['1A', '1B', '1AB', '1AB_ORF', 'S', 'E', 'M', 'N']
    # list of ORFs in final table
    orf_types_final = ['1A', '1B', 'S', 'E', 'M', 'N']


    out_file_name = os.path.splitext(input_file)[0] + '_coord.txt'
    out_file = open(out_file_name, 'w')

    out_file.write(','.join(orf_types_final) + '\n')
    

    dict_coord = {}

    with open(input_file) as handle:
        records = list(SeqIO.parse(handle, 'gb'))
        for rec in records:
            dict_coord[rec.id] = {}
            # counter for polymerase A and B genes
            # if counter==0, haven't met A gene
            pol_count = 0
            print(rec.id)
            for feature in rec.features:
                if feature.type == 'CDS':
                    if 'product' in feature.qualifiers.keys():
                        product = map_feature(feature.qualifiers['product'][0], orf_dict)
                        print(feature.qualifiers['product'][0], product)
                        if product not in dict_coord[rec.id].keys():
                            if product in orf_types_final:
                                dict_coord[rec.id][product] = [int(feature.location._start), int(feature.location._end)]
                            elif product == '1AB':
                                if pol_count == 1:
                                    continue
                                else:
                                    dict_coord[rec.id]['1A'] = [int(feature.location.parts[0]._start), int(feature.location.parts[0]._end)]
                                    dict_coord[rec.id]['1B'] = [int(feature.location.parts[1]._start), int(feature.location.parts[1]._end)]
                                    pol_count += 1
                            elif product == '1AB_ORF':
                                if pol_count == 0:
                                    dict_coord[rec.id]['1A'] = [int(feature.location._start), int(feature.location._end)]
                                    pol_count += 1
                                elif pol_count == 1:
                                    dict_coord[rec.id]['1B'] = [int(feature.location._start), int(feature.location._end)]
                                else:
                                    print('Couldn\'t find annotation \'product\' qualifier for {}'.format(rec.id))
                        
                    elif 'gene' in feature.qualifiers.keys():
                        gene = map_feature(feature.qualifiers['gene'][0], orf_dict)
                        if gene not in dict_coord[rec.id].keys():
                            if gene in orf_types_final:
                                dict_coord[rec.id][gene] = [int(feature.location._start), int(feature.location._end)]
                            elif gene == '1AB':
                                if pol_count == 1:
                                    continue
                                else:
                                    dict_coord[rec.id]['1A'] = [int(feature.location.parts[0]._start), int(feature.location.parts[0]._end)]
                                    dict_coord[rec.id]['1B'] = [int(feature.location.parts[1]._start), int(feature.location.parts[1]._end)]
                                    pol_count += 1
                            elif gene == '1AB_ORF':
                                if pol_count == 0:
                                    dict_coord[rec.id]['1A'] = [int(feature.location._start), int(feature.location._end)]
                                    pol_count += 1
                                elif pol_count == 1:
                                    dict_coord[rec.id]['1B'] = [int(feature.location._start), int(feature.location._end)]
                                else:
                                    print('Couldn\'t find annotation in \'gene\' qualifier for {}'.format(rec.id))
    print(dict_coord)
    out_file.close()


def map_feature(feature, feature_map):
    '''
    feature_map - dictionary, e.g. feature_map['Italia']='ITA'
    feature - possible key from feature_map
    return value if key is in feature_map
    '''
    for k, v in feature_map.items():
        if feature == k:
            return v
    return feature







def retrieve_orf(input_file, coord_file):
    '''
    Splits nucleotide sequences from input_file into 5'NTR,
    one coding region (everything between NTRs), 3'NTR using coordinates of these region from 
    coord_file

    Input:
        input_file - file with nucleotide sequences in fasta-format
        coord_file - text file with coordinates of NTRs
        


    '''
    out_temp = os.path.splitext(input_file)[0]
    
    output_5utr_n = out_temp + '_5utr.fasta'
    output_cds_n = out_temp + '_coding.fasta'
    output_3utr_n = out_temp + '3utr.fasta'
    
    print(output_5utr_n)
    
    list_5utr = list()
    list_cds = list()
    list_3utr = list()

    coord_df = pd.read_csv(coord_file, sep=",", index_col=0)

    seqs = SeqIO.parse(open(input_file), format='fasta')

    for seq in seqs:
        # accession number of sequence
        acc = seq.id.split("_")[0]
        if acc == 'NC' or acc == 'AC':
            acc = '_'.join([acc,seq.id.split("_")[1]])
        if acc in coord_df.index:
            if coord_df['st_5utr'][acc] !=0 and coord_df['e_5utr'][acc]!= 0:
                seq_5utr = copy.deepcopy(seq)
                seq_5utr.seq = seq.seq[:coord_df['e_5utr'][acc]]
                list_5utr.append(seq_5utr)
            seq_cds = copy.deepcopy(seq)
            seq_cds.seq = seq.seq[coord_df['st_cds'][acc]-1:coord_df['e_cds'][acc]]
            list_cds.append(seq_cds)
            
            if coord_df['st_3utr'][acc] !=0 and coord_df['e_3utr'][acc]!= 0:
                seq_3utr = copy.deepcopy(seq)
                seq_3utr.seq = seq.seq[coord_df['st_3utr'][acc]:coord_df['e_3utr'][acc]]
                list_3utr.append(seq_3utr)
        else:
            print('{} not in coord file'.format(acc))
    SeqIO.write(list_5utr, output_5utr_n, "fasta")
    SeqIO.write(list_cds, output_cds_n, "fasta")
    SeqIO.write(list_3utr, output_3utr_n, "fasta")

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-input", "--input_file", type=str,
                        help="Input file", required=True)
    parser.add_argument("-orf_map", "--orf_map_file", type=str,
                        help="Csv-file with short codes for ORFs", required=True)
    args = parser.parse_args()

    orf_coord(args.input_file, args.orf_map_file)