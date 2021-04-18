import argparse
import os
from colour import Color
import seaborn as sns

def color_tree(tree_file_name, color_dict):
    tree_colored_lines = []
    k = 0
    with open(tree_file_name) as tree_file:
        for line in tree_file:
            if line == '\ttaxlabels\n':
                k = 1
                tree_colored_lines.append(line)
                continue
            if line == ';\n' and k == 1:
                k = 0
            if k == 0:
                tree_colored_lines.append(line)
            if k == 1:
                taxa = line.strip('\n').strip('\t')
                already_color = taxa.find('[&!color')
                if already_color != -1:
                    taxa = taxa[:already_color]
                line = '\t' + taxa + '[&!color=' + color_dict[taxa] + ']\n'
                tree_colored_lines.append(line)
    tree_file.close()
    tree_colored = os.path.splitext(tree_file_name)[0] +'_colorgrad.nwk'
    with open(tree_colored, 'w') as tree_colored_file:
        tree_colored_file.writelines(tree_colored_lines)
    tree_colored_file.close()

def create_color_dict(taxa_file_name):
    # reading file with taxa list
    with open(taxa_file_name) as taxa_file:
        list_taxa = [line.strip('\n') for line in taxa_file.readlines()]
    taxa_file.close()

    # create hex codes for gradient colors
    blue = Color("blue")
    green = Color("red")
    grad = list(green.range_to(blue, len(list_taxa)))
    grad_hex = [x.get_hex() for x in grad]
    grad_hex_cor = []
    for hex_code in grad_hex:
        if len(hex_code) == 4:
            hex_code = '#' + 2 * hex_code[1] + 2 * hex_code[2] + 2 * hex_code[3] 
        grad_hex_cor.append(hex_code)
    # create color dict
    dict_color = {}
    for i in range(len(list_taxa)):
        dict_color[list_taxa[i]] = grad_hex_cor[i]
    return dict_color

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-taxa", "--taxa", type=str,
                        help="Text file with taxa names separated by \n", required=True)
    parser.add_argument("-t1", "--tree1", type=str,
                        help="Tree file in nexus format", required=True)
    parser.add_argument("-t2", "--tree2", type=str,
                        help="Tree file in nexus format", required=True)

    args = parser.parse_args()
    # file contains list with original order of taxa
    #tree1 = 'gamma_1-11781.nwk'
    #tree2 = 'gamma_11782-19470.nwk'


    # file contains list with original order of taxa
    #taxa_file_name = 'gamma_1-11781_taxa.txt'

    dict_color = create_color_dict(args.taxa)

    color_tree(args.tree1, dict_color)
    color_tree(args.tree2, dict_color)