import re
import shutil
import csv
import argparse


def csv_reader(csv_file_input, input_str):
    '''
    Reads csv-file and outputs the short name of a variable.

    Input:
        csv_file_input - string - name of csv-file in the directory of the script
        input_str - string - a string value

    Output:
        output_str - string - short designation from csv-file
    '''

    with open(csv_file_input) as csv_file:
        output_str = ''
        reader = csv.DictReader(csv_file, delimiter=",",
                                fieldnames=["base", "new"])
        for line in reader:
            base = re.compile(line["base"])
            if base.match(input_str):
                output_str = line["new"].strip()
        return output_str


def name_tree_leaves(tree, csv):
    renamed_tree = '/'.join(tree.split('/')[:-1]) + '/renamed.tree'
    shutil.copyfile(tree, renamed_tree)
    taxlabels_flag = False
    trees_flag = False
    leaves_list = []
    leaves_list_new = []
    with open(tree, 'r') as in_tree, \
            open(renamed_tree, 'w+') as out_tree:
        for line in in_tree:
            if not (taxlabels_flag or trees_flag) is True:
                out_tree.write(line)
            if re.match('\ttaxlabels', line):
                taxlabels_flag = True
                continue
            if re.match('begin trees;', line):
                trees_flag = True
                continue
            if taxlabels_flag is True:
                if line == ';\n':
                    taxlabels_flag = False
                    trees_flag = False
                    out_tree.write(line)
                    continue
                else:
                    seq_name = re.search(r"[A-Za-z0-9_\-\/.]+", line).group()
                    color = re.search(r'#[0-9a-z]+', line).group()
                    seq_name_new = seq_name + '_' + csv_reader(csv, color)
                    leaves_list.append(seq_name)
                    leaves_list_new.append(seq_name_new)
                    print(seq_name, '--->', seq_name_new)
                if taxlabels_flag is True:
                    out_tree.write(line.split(seq_name)[0] + seq_name_new + line.split(seq_name)[1])
            if trees_flag is True:
                for seq_old, seq_new in zip(leaves_list, leaves_list_new):
                    try:
                        line = (line.split(seq_old)[0] + seq_new + line.split(seq_old)[1])
                    except IndexError:
                        print(seq_old, seq_new)
                        break
                trees_flag = False
                out_tree.write(line)
    print("renamed tree: " + renamed_tree)
    return None


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-tree", "--input_tree", type=str,
                        help="Input tree in nexus format", required=True)
    parser.add_argument("-csv", "--csv_file", type=str,
                        help="Csv-file with colors and corresponding genotypes", required=True)

    args = parser.parse_args()

    name_tree_leaves(args.input_tree, args.csv_file)
