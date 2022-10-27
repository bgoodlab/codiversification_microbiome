# This script calculates the PACo score for the hypothetical example in Fig. 1
# (uses the observed human phylogeny for the E. rectale hosts in the real data)

import sys
import os

permutation = 'all' # 'all', 'country', 'continent'
num_permutations = 1000 # integer
cluster = False # 'True' or 'False'
species = 'Prevotella_copri'

os.system('python create_hypothetical_bacterial_tree.py %s' % species)
os.system('Rscript paco_test.R leylabdata/fake_%s.tre  leylabdata/HostTree.tre leylabdata/HostTree_1.tre leylabdata/HostTree_2.tre %d --permutation=%s --cluster=%s' % (species, num_permutations, permutation, cluster))
	