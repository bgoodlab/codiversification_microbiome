# The purpose of this script is to correct the node names from SNPphylo (which truncates them to 8 characters)

import itertools
from Bio import Phylo
import sys
import numpy

short_long_map = {}
file = open("leylabdata/ID_table_n839.txt","r")
file.readline() # header
for line in file:
	items = line.split()
	tree_id = items[0].strip()
	if len(tree_id)>10:
		short_tree_id = tree_id[:10]
		#print(tree_id, short_tree_id)
	else:
		short_tree_id = tree_id
	
	short_long_map[short_tree_id] = tree_id
		
file.close()

input_filename = sys.argv[1]
output_filename = sys.argv[2]

tree = Phylo.read(input_filename, 'newick')

terminals = tree.get_terminals()
for x in terminals:
	x.name = short_long_map[x.name]
		
Phylo.write(tree, output_filename, 'newick')
