# SaV VP2 tree
python .\get_fragment_alignment.py -input C:\Users\orlov\Documents\term_project\sapovirus_genomes_orf2.fasta -ref C:\Users\orlov\Documents\term_project\SaV_reference_2.fasta -frag_s 6851 -frag_e 7349 -output C:\Users\orlov\Desktop\tmp\sapovirus_vp2_aln.fasta

python ..\sapovirus\genotyping.py -in_rep C:\Users\orlov\Desktop\tmp\ -in_tree C:\Users\orlov\Desktop\tmp\CP_tree_bootstrap1000_colored_nodes.nwc -in_csv C:\Users\orlov\Desktop\tmp\Color_codes.txt

iqtree -s C:\Users\orlov\Desktop\tmp\sapovirus_vp2_aln_genotyped.fasta -bb 1000