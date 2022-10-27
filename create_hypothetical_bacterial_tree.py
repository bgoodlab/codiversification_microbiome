#
# The purpose of this script is to create a fake bacterial tree corresponding
# to the hypothetical example shown in Fig. 1.
#

import pylab
import itertools
from Bio import Phylo
from Bio.Phylo import TreeConstruction
import sys
import numpy
import statsmodels.api as sm

import matplotlib.colors as colors
import matplotlib.cm as cmx
from math import log10,ceil
import matplotlib as mpl
import matplotlib.gridspec as gridspec
from numpy.random import randint,random_sample
import matplotlib.colors as mcolors

sample_country_map = {}
sample_continent_map = {}
allowed_geographic_samples = []
file = open("leylabdata/ID_table_n839.txt","r")
file.readline() # header
for line in file:
	items = line.split()
	#print(items)
	sample = items[0].strip()
	country = items[5].strip()
	continent = items[6].strip()
	
	sample_country_map[sample] = country
	sample_continent_map[sample] = continent
	
	if country=='Gabon': # or continent=='Africa':
		allowed_geographic_samples.append(sample)
	#if continent=='Europe':
	#	allowed_geographic_samples.append(sample)	
file.close()

#allowed_geographic_samples = set(allowed_geographic_samples)
allowed_geographic_samples = set(sample_country_map.keys())

species = sys.argv[1]
bacteria_filename = "leylabdata/Adult_RAxML_bipartitions.s__%s.StrainPhlAn3.tre" % species

fake_filename = "leylabdata/fake_%s.tre" % species


bacteria_tree = Phylo.read(bacteria_filename, 'newick') 
bacteria_terminals = bacteria_tree.get_terminals()

strain_probability_map = {'Africa':0.7,'Europe':0.6,'Asia':0.3}
# test with different probabilities to ensure that it is not that sensitive to the first choice. 
#strain_probability_map = {'Africa':0.7,'Europe':0.3,'Asia':0.6}

do = 2.0
di = 0.0

bacteria_names = []
bacteria_continents = []
bacteria_strains = []

alpha_bacteria = []
beta_bacteria = []

for x in bacteria_terminals:
	name = x.name
	bacteria_names.append(name)
	continent = sample_continent_map[name]
	bacteria_continents.append(continent)
	
	p = random_sample()
	strain = 0
	if p<strain_probability_map[continent]:
		strain = 1
	bacteria_strains.append(strain)
	
	if strain==1:
		alpha_bacteria.append(name)
	else:
		beta_bacteria.append(name)

alpha_string = ",".join(["%s:%g" % (name,di) for name in alpha_bacteria])
beta_string = ",".join(["%s:%g" % (name,di) for name in beta_bacteria])
total_string = "((%s):%g,(%s):%g);" % (alpha_string, do, beta_string, do)

print(len(bacteria_names),len(alpha_bacteria),len(beta_bacteria))
file = open(fake_filename,"w")
file.write(total_string)
file.close()
