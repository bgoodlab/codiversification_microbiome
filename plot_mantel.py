import pylab
import itertools
from Bio import Phylo
import sys
import numpy
import statsmodels.api as sm
import matplotlib.colors as mcolors
from numpy.random import shuffle
from scipy.stats import linregress

sample_country_map = {}
sample_continent_map = {}
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
		
file.close()

bacteria_filename = sys.argv[1]
host_filename = sys.argv[2]
species = bacteria_filename
species2 = host_filename

bacteria_tree = Phylo.read(bacteria_filename, 'newick')
host_tree = Phylo.read(host_filename,'newick')

bacteria_terminals = bacteria_tree.get_terminals()
host_terminals = host_tree.get_terminals()

bacteria_names=[]
for x in bacteria_terminals:
	bacteria_names.append(x.name)
bacteria_names = set(bacteria_names)

host_names=[]
for x in host_terminals:
	host_names.append(x.name)
host_names = set(host_names)

allowed_names = (bacteria_names & host_names)

new_allowed_names = set()
for sample in allowed_names:
	if (sample_country_map[sample]=='Gabon'):
		new_allowed_names.add(sample)

# Restrict to certain subset or not
#allowed_names = new_allowed_names

# get distance dictionary of bacteria
print("Loading bacteria distances from %d leaves..." % len(allowed_names))
bacteria_distance_dict = {}
bacteria_names = []
num_processed = 0
for x in bacteria_terminals:
	if x.name not in allowed_names:
		continue
	num_processed+=1
	if num_processed % 100 == 0:
		print("Processed %d/%d" % (num_processed, len(bacteria_terminals)))
	bacteria_names.append(x.name)
	bacteria_distance_dict[x.name]={}
	for y in bacteria_terminals:
		if y.name not in allowed_names:
			continue
		if x==y:
			d=0
		else:
			d = bacteria_tree.distance(x, y)
		bacteria_distance_dict[x.name][y.name] = d
print("Done!")	  

# get distance dictionary of host
print("Loading host distances...")
host_distance_dict = {}
for x in host_terminals:
	if x.name not in allowed_names:
		continue
	host_distance_dict[x.name]={}
	for y in host_terminals:
		if y.name not in allowed_names:
			continue
		if x==y:
			d=0
		else:
			d = host_tree.distance(x, y)
		host_distance_dict[x.name][y.name] = d
print("Done!")

# Now turn it into distance matrix
bacteria_ds = []
host_ds = []
for i in range(0,len(bacteria_names)):
	bacteria_row = []
	host_row = []
	for j in range(0,len(bacteria_names)):
		host_row.append(host_distance_dict[bacteria_names[i]][bacteria_names[j]])
		bacteria_row.append(bacteria_distance_dict[bacteria_names[i]][bacteria_names[j]])
	
	bacteria_ds.append(bacteria_row)
	host_ds.append(host_row)
	
host_ds = numpy.array(host_ds)
bacteria_ds = numpy.array(bacteria_ds)

#smoothed = sm.nonparametric.lowess(exog=host_ds, endog=bacteria_ds, frac=0.2)

pylab.plot(host_ds,bacteria_ds,'k.',alpha=0.1,markersize=2)
#pylab.hist2d(host_ds, bacteria_ds, bins=(50, 50),norm=mcolors.PowerNorm(0.3),cmin=1,cmap="copper_r")

def mantel_test(D1,D2,num_bootstraps=1000):
	
	sample_idxs = numpy.arange(0,D1.shape[0])
	
	triu_indices = numpy.triu_indices(D1.shape[0],1)
	
	bootstrapped_D1 = numpy.array(D1,copy=True)

	LR = linregress(D1[triu_indices],D2[triu_indices])
	m = LR.slope
	observed_correlation = LR.rvalue 
	
	#observed_correlation = numpy.corrcoef(D1[triu_indices],D2[triu_indices])[0,1]
	bootstrapped_correlations = []
	for b in range(0,num_bootstraps):
		shuffle(sample_idxs)
		bootstrapped_D1[:,:] = bootstrapped_D1[sample_idxs,:]
		bootstrapped_D1[:,:] = bootstrapped_D1[:,sample_idxs]
		
		#print("Should be zero:", numpy.diag(bootstrapped_D1))
		
		bootstrapped_correlation = numpy.corrcoef(bootstrapped_D1[triu_indices],D2[triu_indices])[0,1]
		bootstrapped_correlations.append(bootstrapped_correlation)
		
	bootstrapped_correlations = numpy.array(bootstrapped_correlations)
	
	pvalue = ((bootstrapped_correlations>=observed_correlation).sum()+1.0)/(len(bootstrapped_correlations)+1.0)
	
	return m, observed_correlation, pvalue
		
	
m, observed_correlation, pvalue = mantel_test(host_ds, bacteria_ds)

print(m, observed_correlation, pvalue)

pylab.savefig('mantel_correlation.png',bbox_inches='tight',dpi=300)
