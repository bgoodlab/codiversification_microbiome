import pylab
import numpy
import sys
import os

permutation = sys.argv[1] # 'all', 'country', 'continent'
num_permutations = int(sys.argv[2]) # integer
cluster = sys.argv[3] # 'True' or 'False'
speciess = []

file = open("paco_results.csv","r")
for line in file:
    items = line.split(",")
    species = items[0].strip()
    es = float(items[1])
    qvalue = float(items[2])

    speciess.append(species)

if cluster=='True':
	output_filename = "paco_output_%s_%d_clusterhumans.txt" % (permutation, num_permutations)
else:
	output_filename = "paco_output_%s_%d_nocluster.txt" % (permutation, num_permutations)
	
os.system('touch %s' % output_filename)
os.system('echo > %s' % output_filename)
for i in range(0,len(speciess)):
	species = speciess[i]
	print("Running %d/%d (%s)..." % (i+1,len(speciess),species))
	os.system('Rscript leylabcode/ben_paco.R leylabdata/Adult_RAxML_bipartitions.s__%s.StrainPhlAn3.tre leylabdata/HostTree.tre leylabdata/HostTree_1.tre leylabdata/HostTree_2.tre %d --permutation=%s --cluster=%s >> %s 2>/dev/null' % (species, num_permutations, permutation, cluster, output_filename))
	