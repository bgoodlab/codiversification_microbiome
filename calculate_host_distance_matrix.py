import pylab
import itertools
from Bio import Phylo
import sys
import numpy
import statsmodels.api as sm


# Load host genotype matrix
host_filename = sys.argv[1] 
host_genotype_matrix = []
file = open(host_filename,"r")
samples = [item.strip() for item in file.readline().split(",")]
for line in file:
	items = line.split(",")
	genotypes = numpy.array([float(item) for item in items])
	
	# Polarize to minor allele
	k = genotypes[genotypes>-0.5].sum()
	n = (genotypes>-0.5).sum()
	if k > n/2:
		genotypes = 1-genotypes
		genotypes[genotypes>1.5] = -1
	
	host_genotype_matrix.append(genotypes)
genotype_matrix = numpy.array(host_genotype_matrix)

good_idxs = (genotype_matrix>-0.5)
ns = good_idxs.sum(axis=1)
ks = (genotype_matrix*good_idxs).sum(axis=1)
fs = ks*1.0/ns

joint_good_idxs = good_idxs[:,:,None]*good_idxs[:,None,:]
	
scaled_genotype_matrix = (genotype_matrix-fs[:,None])/numpy.sqrt((fs*(1-fs))[:,None])
	
D = (scaled_genotype_matrix[:,:,None]*scaled_genotype_matrix[:,None,:]*joint_good_idxs).sum(axis=0)
D = D / joint_good_idxs.sum(axis=0)

print(", ".join(samples))
for i in range(0,D.shape[0]):
	print(", ".join([str(item) for item in D[i,:]]))
