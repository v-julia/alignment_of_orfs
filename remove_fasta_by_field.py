from Bio import SeqIO
import re
import argparse


def remove_fasta_by_field(in_f, field, out_f):
    with open(in_f, "r") as handle,\
         open(out_f, "w") as output:
        if "/" in field:
            field = field.replace("/", "\\")
        for rec in SeqIO.parse(handle, "fasta"):
            if not re.fullmatch(rf"{field}", rec.id):
                SeqIO.write(rec, output, "fasta")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-input", "--input_file", type=str,
                        help="Input alignment in fasta", required=True)
    parser.add_argument("-field", "--field", type=str,
                        help="The field by which fasta will be removed (regexp)")
    parser.add_argument("-output", "--output_file", type=str,
                        help="Output alignment in fasta")
    args = parser.parse_args()
remove_fasta_by_field(args.input_file, args.field, args.output_file)
