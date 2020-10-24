import argparse
import os

from Bio.Seq import Seq
from Bio import SeqIO
from Bio import AlignIO

def cut_fragment(ref, start, end):
    rec_ref_name = os.path.splitext(ref)[0] + '_fragment.fasta'
    rec_ref = next(SeqIO.parse(open(ref), 'fasta'))
    seq_frag = (rec_ref).seq[start:end]
    setattr(rec_ref, 'seq', seq_frag)
    SeqIO.write(rec_ref, rec_ref_name, 'fasta')
    return rec_ref
    
def aln_sample(input_file, rec_ref):
    sample_name = os.path.splitext(input_file)[0] + '_frag_sample.fasta'
    sample = SeqIO.parse(open(input_file), 'fasta')
    sample_list = []
    sample_list.append(rec_ref)
    for record in sample:
        print(record)
        record.seq = Seq((str(record.seq)).replace('-', ''))
        sample_list.append(record)
        print(record)
    SeqIO.write(sample_list, sample_name, 'fasta')
    alignment_command = 'python .\\trans_alignment.py -input {}'.format(sample_name)

    os.system(alignment_command)
    aln_name = os.path.splitext(sample_name)[0] + '_tr_aln_rt.fasta'
    return aln_name

def cut_aln(aln_name, output_file):
    aln = AlignIO.read(open(aln_name), 'fasta')
    ref_seq = aln[0].seq
    count_start = -1
    for pos in ref_seq:
        count_start += 1
        if pos != '-':
            start_pos = count_start
            break
    count_end = len(ref_seq) + 1
    for pos in reversed(ref_seq):
        count_end -= 1
        if pos != '-':
            end_pos = count_end
            break
    frag_aln_list = []
    for rec in SeqIO.parse(open(aln_name), 'fasta'):
        setattr(rec, 'seq', rec.seq[start_pos:end_pos])
        frag_aln_list.append(rec)
    del(frag_aln_list[0]) #remove reference
    SeqIO.write(frag_aln_list, output_file, 'fasta')
    return
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-input", "--input_file", type=str,
                        help="Input file in fasta format", required=True)
    parser.add_argument("-ref", "--reference", type=str,
                        help="Reference (one record) in fasta format", required=True)
    parser.add_argument("-frag_s", "--fragment_start", type=int,
                        help="The beginning of the required fragment in reference sequence ", required=True)
    parser.add_argument("-frag_e", "--fragment_end", type=int,
                        help="The end of the required fragment in reference sequence ", required=True)
    parser.add_argument("-output", "--output_file", type=str,
                        help="Output file in fasta format", required=True)

    args = parser.parse_args()

    rec_ref = cut_fragment(args.reference, args.fragment_start, args.fragment_end)
    aln_name = aln_sample(args.input_file, rec_ref)
    cut_aln(aln_name, args.output_file)