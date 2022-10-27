import pylab
import itertools
from Bio import Phylo
import sys
import numpy
import statsmodels.api as sm




if len(sys.argv) < 2:
    species = "host"
else:
    species = sys.argv[1]

# Load host genotype matrix
host_filename = "host_distance_matrix_all.txt"
host_distance_matrix = []
file = open(host_filename,"r")
samples = [item.strip() for item in file.readline().split(",")]
for line in file:
	items = line.split(",")
	host_distance_matrix.append( numpy.array([1-float(item) for item in items]) )
	
host_distance_matrix = numpy.array(host_distance_matrix)
allowed_names = set(samples)
host_filename = "leylabdata/HostTree.tre"
host_tree = Phylo.read(host_filename,'newick')
host_terminals = host_tree.get_terminals()

tree_names = set([x.name for x in host_terminals])

# Check for perfect match
print(len(allowed_names & tree_names), len(tree_names))
print(allowed_names - tree_names)
print(tree_names - allowed_names)

print("Loading tree...")
num_processed=0
host_distance_dict={}
for i in range(0,len(host_terminals)):
	x = host_terminals[i]
	if x.name not in allowed_names:
		print("Shouldn't happen!",x.name)
		continue
	num_processed+=1
	if num_processed % 100 == 0:
		print("Processed %d/%d" % (num_processed, len(host_terminals)))
    
	host_distance_dict[x.name]={}
	for j in range(i,len(host_terminals)):
		y = host_terminals[j]
    	
		if y.name not in allowed_names:
			print("Shouldn't happen!", x.name, y.name)
			continue
		if x==y:
			d=0
		else:
			d = host_tree.distance(x, y)
            
		host_distance_dict[(x.name, y.name)] = d
		host_distance_dict[(y.name, x.name)] = d
		

print("Done!")



# Now turn it into distance list
pca_ds = []
tree_ds = []
for i in range(0,len(samples)):
    for j in range(i+1,len(samples)):
    	if (samples[i],samples[j]) not in host_distance_dict:
    		print("j not in dict", i,j, samples[i],samples[j])
    	tree_ds.append(host_distance_dict[(samples[i],samples[j])])
    	pca_ds.append(host_distance_matrix[i,j])

pylab.plot(tree_ds,pca_ds,'k.',alpha=0.5,markersize=2)
pylab.xlabel("Host tree distance")
pylab.ylabel("Host PCA distance")
pylab.savefig('host_pca_tree_correlation.png',bbox_inches='tight',dpi=300)
