import itertools
from Bio import Phylo
import sys
import numpy


species = sys.argv[1]

bacteria_filename = "leylabdata/Adult_RAxML_bipartitions.s__%s.StrainPhlAn3.tre" % species
host_filename = "leylabdata/HostTree.tre"
new_host_filename = "leylabdata/HostTree_%s.tre" % species

bacteria_tree = Phylo.read(bacteria_filename, 'newick')
host_tree = Phylo.read(host_filename,'newick')

bacteria_terminals = bacteria_tree.get_terminals()
host_terminals = host_tree.get_terminals()

bacteria_names=[]
for x in bacteria_terminals:
    bacteria_names.append(x.name)

host_subtree = host_tree.common_ancestor(bacteria_names)

output_file = open(new_host_filename, 'w')
Phylo.write(host_subtree, output_file, "newick")
output_file.close()
