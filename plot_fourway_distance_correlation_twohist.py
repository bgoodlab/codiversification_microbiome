import pylab
import itertools
from Bio import Phylo
import sys
import numpy
import statsmodels.api as sm

import matplotlib.colors as colors
import matplotlib.cm as cmx
from math import log10,ceil
import matplotlib as mpl
import matplotlib.gridspec as gridspec
from numpy.random import randint
import matplotlib.colors as mcolors

other_color = '#3B75AF'
continent_color = '#519E3E'
country_color = '#EF8636'
related_color = '0.7'

def get_pretty_name(species_name):
	items = species_name.split("_")
	
	genus = items[0]
	species = items[1]
	
	if species=='thetaiotaomicron':
		species='theta' # usually abbreviated as B. theta
	return "%s. %s" % (genus[0],species)

mpl.rcParams['font.size'] = 7
mpl.rcParams['lines.linewidth'] = 0.5
mpl.rcParams['legend.frameon']  = False
mpl.rcParams['legend.fontsize']  = 'small'

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

species = {}
species[1] = sys.argv[1]
species[2] = sys.argv[2]

bacteria_filename = {}
bacteria_filename[1] = "leylabdata/Adult_RAxML_bipartitions.s__%s.StrainPhlAn3.tre" % species[1]
bacteria_filename[2] = "leylabdata/Adult_RAxML_bipartitions.s__%s.StrainPhlAn3.tre" % species[2]

host_filename = {}
host_filename[0] = "leylabdata/HostTree.tre"
host_filename[1] = "leylabdata/HostTree_1.tre"
host_filename[2] = "leylabdata/HostTree_2.tre"

bacteria_tree = {s: Phylo.read(bacteria_filename[s], 'newick') for s in [1,2]}
host_tree = {h: Phylo.read(host_filename[h],'newick') for h in [0,1,2]}

bacteria_terminals = {s: bacteria_tree[s].get_terminals() for s in [1,2]}
host_terminals = {h: host_tree[h].get_terminals() for h in [0,1,2]}

bacteria_names = {}
for i in [1,2]:
	bacteria_names[i] = []
	for x in bacteria_terminals[i]:
		bacteria_names[i].append(x.name)
	bacteria_names[i] = set(bacteria_names[i])

host_names = {}
for i in [0,1,2]:
	host_names[i]=[]
	for x in host_terminals[i]:
		host_names[i].append(x.name)
	host_names[i] = set(host_names[i])

allowed_names = {}
allowed_names[1] = (bacteria_names[1] & host_names[0] & host_names[1] & host_names[2]) & allowed_geographic_samples
allowed_names[2] = (bacteria_names[2] & host_names[0] & host_names[1] & host_names[2]) & allowed_geographic_samples
allowed_names_or = ((bacteria_names[1] | bacteria_names[2]) & host_names[0] & host_names[1] & host_names[2]) & allowed_geographic_samples
allowed_names_and = ((bacteria_names[1] & bacteria_names[2]) & host_names[0] & host_names[1] & host_names[2]) & allowed_geographic_samples

# get distance dictionary of bacteria
bacteria_distance_dict = {}
for s in [1,2]:
	print("Loading %s..." % species[s])
	bacteria_distance_dict[s]={}
	num_processed = 0
	for i in range(0,len(bacteria_terminals[s])):
		x = bacteria_terminals[s][i]
		if x.name not in allowed_names[s]:
			continue
		num_processed+=1
		if num_processed % 100 == 0:
			print("Processed %d/%d" % (num_processed, len(allowed_names[s])))
		
		for j in range(i,len(bacteria_terminals[s])):
			y = bacteria_terminals[s][j]
			if y.name not in allowed_names[s]:
				continue
			if x==y:
				d=0
			else:
				d = bacteria_tree[s].distance(x, y)
			
			bacteria_distance_dict[s][(x.name,y.name)] = d
			bacteria_distance_dict[s][(y.name,x.name)] = d
		
	print("Done!")	  

host_distance_dict = {}
for h in [0,1,2]:
	# get distance dictionary of host
	print("Loading host %d distances..." % h)
	host_distance_dict[h] = {}
	for i in range(0,len(host_terminals[h])):
		x = host_terminals[h][i]
		if x.name not in allowed_names_or:
			continue
		for j in range(i,len(host_terminals[h])):
			y = host_terminals[h][j] 
			if y.name not in allowed_names_or:
				continue
			if x==y:
				d=0
			else:
				d = host_tree[h].distance(x, y)
			host_distance_dict[h][(x.name,y.name)] = d
			host_distance_dict[h][(y.name,x.name)] = d
	print("Done!")

# Now turn it into distance list
bacteria_ds = {}
host_ds = {}
host_host_ds = {h:[] for h in [1,2]}
bacteria_dmatrices = {}
host_dmatrices = {}

same_country = {}
same_continent = {}
same_continent_diff_country = {}
diff_continent = {}

host_same_country = []
host_same_continent = []

host_host_samples = set()
for s in [1,2]:
	samples = list(allowed_names[s])
	bacteria_ds[s] = []
	host_ds[s] = []
	same_country[s] = []
	same_continent[s] = []
	bacteria_dmatrices[s] = numpy.zeros((len(samples),len(samples)))*1.0
	host_dmatrices[s] = numpy.zeros((len(samples),len(samples)))*1.0
	
	for i in range(0,len(samples)):
		for j in range(i+1,len(samples)):
			
			bacteria_d = bacteria_distance_dict[s][(samples[i],samples[j])]
			host_d = host_distance_dict[0][(samples[i],samples[j])] # Use full host phylogenetic tree to maximize correlation between bacteria and host
			
			bacteria_ds[s].append(bacteria_d)
			host_ds[s].append(host_d)
			
			bacteria_dmatrices[s][i,j] = bacteria_d
			bacteria_dmatrices[s][j,i] = bacteria_d
			
			host_dmatrices[s][i,j] = host_d
			host_dmatrices[s][j,i] = host_d
			
			is_same_country = sample_country_map[samples[i]]==sample_country_map[samples[j]]
			is_same_continent = sample_continent_map[samples[i]]==sample_continent_map[samples[j]]
			
			same_country[s].append(is_same_country)
			same_continent[s].append(is_same_continent)
			
			if (samples[i],samples[j]) not in host_host_samples:
				for h in [1,2]:
					host_host_ds[h].append(host_distance_dict[h][(samples[i],samples[j])])
				host_same_country.append(is_same_country)
				host_same_continent.append(is_same_continent)
				host_host_samples.add((samples[i],samples[j]))
				host_host_samples.add((samples[j],samples[i]))
			
	bacteria_ds[s] = numpy.array(bacteria_ds[s])
	host_ds[s] = numpy.array(host_ds[s])
	same_country[s] = numpy.array(same_country[s])
	same_continent[s] = numpy.array(same_continent[s])
	same_continent_diff_country[s] = numpy.logical_and(same_continent[s], numpy.logical_not(same_country[s]))
	diff_continent[s] = numpy.logical_not(same_continent[s])

for h in [1,2]:	
	host_host_ds[h] = numpy.array(host_host_ds[h])
host_same_country = numpy.array(host_same_country)
host_same_continent = numpy.array(host_same_continent)
host_same_continent_diff_country = numpy.logical_and(host_same_continent, numpy.logical_not(host_same_country))
host_diff_continent = numpy.logical_not(host_same_continent)

# Now do host vs host distance matrix
host_host_dmatrices = {}
samples = list(allowed_names_or)
for h in [1,2]:
	host_host_dmatrices[h] = numpy.zeros((len(samples),len(samples)))*1.0
	for i in range(0,len(samples)):
		for j in range(i+1,len(samples)):
			d = host_distance_dict[h][(samples[i],samples[j])]
			host_host_dmatrices[h][i,j] = d
			host_host_dmatrices[h][j,i] = d
				

# Now get merged distance list for just the two bacteria
bac_bac_ds = {s:[] for s in [0,1,2]}
samples = list(allowed_names_and)
bac_bac_dmatrices = {s:numpy.zeros((len(samples),len(samples))) for s in [0,1,2]}
bac_bac_same_country = []
bac_bac_same_continent = []
for i in range(0,len(samples)):
	for j in range(i+1,len(samples)):
		for s in [1,2]:
			d = bacteria_distance_dict[s][(samples[i],samples[j])]
			bac_bac_ds[s].append(d)
			bac_bac_dmatrices[s][i,j] = d
			bac_bac_dmatrices[s][j,i] = d
		
		bac_bac_dmatrices[0][i,j] = host_distance_dict[s][(samples[i],samples[j])]
		bac_bac_dmatrices[0][j,i] = host_distance_dict[s][(samples[i],samples[j])]
		
		bac_bac_ds[0].append( host_distance_dict[s][(samples[i],samples[j])] )	
		is_same_country = sample_country_map[samples[i]]==sample_country_map[samples[j]]
		is_same_continent = sample_continent_map[samples[i]]==sample_continent_map[samples[j]]
			
		bac_bac_same_country.append(is_same_country)
		bac_bac_same_continent.append(is_same_continent)

bac_bac_same_country = numpy.array(bac_bac_same_country)
bac_bac_same_continent = numpy.array(bac_bac_same_continent)
bac_bac_same_continent_diff_country = numpy.logical_and(bac_bac_same_continent, numpy.logical_not(bac_bac_same_country))
bac_bac_diff_continent = numpy.logical_not(bac_bac_same_continent)

for s in [0,1,2]:
	bac_bac_ds[s] = numpy.array(bac_bac_ds[s])
	

max_host_d = max([host_host_ds[h].max() for h in [1,2]]) 
max_bacteria_d = {s: bacteria_ds[s].max() for s in [1,2]}

fig = pylab.figure(figsize=(7, 6))

# Set up grids to hold figure panels
outer_grid = gridspec.GridSpec(2, 2, height_ratios=[1,1], width_ratios = [1,1], hspace=0.3, wspace=0.3)

bacteria_axis = {}
bacteria_hist_axis = {}
bacteria_host_hist_axis = {}

for idx in range(1,4):
	
	if idx==3:
		grid_loc = outer_grid[0,1]
	elif idx==2:
		grid_loc = outer_grid[1,1]
	else:
		grid_loc = outer_grid[0,0]
		
	inner_grid = gridspec.GridSpecFromSubplotSpec(2, 2,
                width_ratios=[1,0.25],
                height_ratios=[0.25,1],
                wspace=0.05,
                hspace=0.05,
                subplot_spec=grid_loc)
                
	if idx==3:
		host_axis = pylab.Subplot(fig, inner_grid[1,0])
		fig.add_subplot(host_axis)
		host_hist_axis = pylab.Subplot(fig, inner_grid[1,1])
		fig.add_subplot(host_hist_axis)
		host_host_hist_axis = pylab.Subplot(fig, inner_grid[0,0])
		fig.add_subplot(host_host_hist_axis)

	else:

		bacteria_axis[idx] = pylab.Subplot(fig, inner_grid[1,0])
		fig.add_subplot(bacteria_axis[idx])
		
		bacteria_hist_axis[idx] = pylab.Subplot(fig, inner_grid[1,1])
		fig.add_subplot(bacteria_hist_axis[idx])
		
		bacteria_host_hist_axis[idx] = pylab.Subplot(fig, inner_grid[0,0])
		fig.add_subplot(bacteria_host_hist_axis[idx])
		
inner_grid = gridspec.GridSpecFromSubplotSpec(3, 3,
                width_ratios=[0.1,1,0.15],
                height_ratios=[0.1,1,0.1],
                wspace=0.05,
                hspace=0.05,
                subplot_spec=outer_grid[1,0])		
temp_axis = pylab.Subplot(fig, inner_grid[1,1])
fig.add_subplot(temp_axis)

host_axis.set_ylabel('Host phylogenetic distance (2/2)')
host_axis.set_xlabel('Host phylogenetic distance (1/2)')


fig2 = pylab.figure(figsize=(3.42,2.5))

# Set up grids to hold figure panels
outer_grid = gridspec.GridSpec(2, 2,width_ratios=[1,0.25],height_ratios=[0.25,1],wspace=0.05,hspace=0.05)
bac_bac_axis = pylab.Subplot(fig2, outer_grid[1,0])
fig2.add_subplot(bac_bac_axis)
bac_bac_hist_axis = pylab.Subplot(fig2, outer_grid[1,1])
fig2.add_subplot(bac_bac_hist_axis)
bac_bac_bac_hist_axis = pylab.Subplot(fig2, outer_grid[0,0])
fig2.add_subplot(bac_bac_bac_hist_axis)

fig3 = pylab.figure(figsize=(3.42,3))

# Set up grids to hold figure panels
outer_grid = gridspec.GridSpec(1, 1)
other_bac_bac_axis = pylab.Subplot(fig3, outer_grid[0])
fig3.add_subplot(other_bac_bac_axis)


	
for s in reversed(range(1,5)):
	print("Procesing index %d" % s)
	if s==3: # the host host one
		axis = host_axis
		hist_axis = host_hist_axis
		hhist_axis = host_host_hist_axis
		dxs = host_host_ds[1]
		dys = host_host_ds[2]
		country_idxs = host_same_country
		continent_idxs = host_same_continent_diff_country
		other_idxs = host_diff_continent
		max_d = max_host_d
		max_dx = max_host_d
		max_dy = max_host_d
		
		dmatrixx = host_host_dmatrices[1]
		dmatrixy = host_host_dmatrices[2]
		
	elif s==4: # the bac bac one
		axis = bac_bac_axis
		hist_axis = bac_bac_hist_axis
		hhist_axis = bac_bac_bac_hist_axis
		dys = bac_bac_ds[1]
		dxs = bac_bac_ds[2]
		country_idxs = bac_bac_same_country
		continent_idxs = bac_bac_same_continent_diff_country
		other_idxs = bac_bac_diff_continent
		max_d = max_bacteria_d[2]
		
		max_dx = max_bacteria_d[2]
		max_dy = max_bacteria_d[1]
		
		dmatrixy = bac_bac_dmatrices[1]
		dmatrixx = bac_bac_dmatrices[2]
		dmatrixz = bac_bac_dmatrices[0]
		
	else:
		axis = bacteria_axis[s]
		hist_axis = bacteria_hist_axis[s]
		hhist_axis = bacteria_host_hist_axis[s]
		dxs = host_ds[s]
		dys = bacteria_ds[s]
		country_idxs = same_country[s]
		continent_idxs = same_continent_diff_country[s]
		other_idxs = diff_continent[s]
		max_d = max_bacteria_d[s]
		
		max_dx = max_host_d
		max_dy = max_bacteria_d[s]
		
		dmatrixx = host_dmatrices[s]
		dmatrixy = bacteria_dmatrices[s]
		
	#line1, = axis.plot(dxs[other_idxs], dys[other_idxs], '.',alpha=0.25,markersize=1)
	#line2, = axis.plot(dxs[continent_idxs], dys[continent_idxs], '.',alpha=0.25,markersize=1)
	#line3, = axis.plot(dxs[country_idxs], dys[country_idxs], '.',alpha=0.25,markersize=1)
		
	axis.hist2d(dxs, dys, bins=(50, 50),norm=mcolors.PowerNorm(0.3),cmin=1,cmap="copper_r")
	
	line3, = axis.plot([-2],[-2],'o',label='Diff continent',markersize=3,color=other_color)
	line2, = axis.plot([-2],[-2],'o',label='Diff country',markersize=3,color=continent_color)
	line1, = axis.plot([-2],[-2],'o',label='Same country',markersize=3,color=country_color)
	line0, = axis.plot([-2],[-2],'s',label='Most related',markersize=3,color=related_color)
	
		
	hist_axis.hist(dys[other_idxs], bins=(numpy.linspace(0,max_dy,50)-0.001),histtype='step',density=True,orientation='horizontal',color=other_color)
	hist_axis.hist(dys[continent_idxs], bins=(numpy.linspace(0,max_dy,50)+0.001), histtype='step',density=True,orientation='horizontal',color=continent_color)
	hist_axis.hist(dys[country_idxs], bins=numpy.linspace(0,max_dy,50),histtype='step',density=True,orientation='horizontal',color=country_color)
	
	if s==4:
		other_bac_bac_axis.scatter(dxs[country_idxs], dys[country_idxs],c=bac_bac_ds[0][country_idxs],alpha=0.25,cmap='plasma',marker='.',s=1)
		
		most_related = []
		for i in range(0,dmatrixy.shape[0]):
			ds = numpy.sort(dmatrixz[i,:])
			#idx = min([20,len(dmatrixy.shape[0])])
			dmax = ds[1]
			good_idxs = (dmatrixz[i,:]>0)*(dmatrixz[i,:]<=dmax)*(dmatrixz[i,:]<=0.4) 
			most_related.extend(dmatrixy[i,good_idxs])
	
		#hist_axis.hist(most_related, bins=numpy.linspace(0,1,50),histtype='step',density=True,orientation='horizontal')
	
	
	# Now plot only those most related (up to ~0.4, which is the top amount within country)
	most_related = []
	for i in range(0,dmatrixy.shape[0]):
		ds = numpy.sort(dmatrixx[i,:])
		#idx = min([20,len(dmatrixy.shape[0])])
		dmax = ds[1]
		good_idxs = (dmatrixx[i,:]>0)*(dmatrixx[i,:]<=dmax)*(dmatrixx[i,:]<=0.4) 
		most_related.extend(dmatrixy[i,good_idxs])
	
	hist_axis.hist(most_related, bins=numpy.linspace(0,max_dy,50),histtype='stepfilled',density=True,orientation='horizontal',alpha=0.5,color=related_color)
	
	# Now do same thing but for bacteria strain
	hhist_axis.hist(dxs[other_idxs], bins=(numpy.linspace(0,max_dx,50)-0.001),histtype='step',density=True,color=other_color)
	hhist_axis.hist(dxs[continent_idxs], bins=(numpy.linspace(0,max_dx,50)+0.001), histtype='step',density=True,color=continent_color)
	hhist_axis.hist(dxs[country_idxs], bins=numpy.linspace(0,max_dx,50),histtype='step',density=True,color=country_color)
	
	
	
	most_related = []
	for i in range(0,dmatrixy.shape[0]):
		ds = numpy.sort(dmatrixy[i,:])
		#idx = min([20,len(dmatrixy.shape[0])])
		dmax = ds[1]
		good_idxs = (dmatrixy[i,:]>0)*(dmatrixy[i,:]<=dmax)
		most_related.extend(dmatrixx[i,good_idxs])
	
	hhist_axis.hist(most_related, bins=numpy.linspace(0,max_dx,50),histtype='stepfilled',density=True,alpha=0.5,color=related_color)
	
	# Now do same thing but for strains as closely related as within country
	good_idxs = (dys<=numpy.median(dys[country_idxs]))
	#hhist_axis.hist(dxs[good_idxs], bins=numpy.linspace(0,1,50),histtype='step',density=True)
	
	if s==3:
		axis.legend(loc='upper left',frameon=False,numpoints=1,handletextpad=0.1)

# Now plot temperature correlation
file = open("phenotypic_correlations.csv")
speciess = []
ess = []
lowtemp_growths = []
qvalues = []
genera = []
phyla = []
file.readline() # header
for line in file:
	items = line.split(",")
	species_name = items[0].strip()
	es = float(items[4])
	phylum = items[3].strip()
	lowtemp_growth = float(items[6])
	qvalue = float(items[7])
	
	#genus = (species.split("_")[0]).split()[0]
	genus = items[8].strip() # actually family
	speciess.append(species_name)
	ess.append(es)
	lowtemp_growths.append(lowtemp_growth)
	phyla.append(phylum)
	genera.append(genus)
	qvalues.append(qvalue)
	 

allowed_genera = set(genera)

for genus in sorted(allowed_genera):
	xs = []
	ys = []
	subqvalues = []
	subspeciess = []
	for i in range(0,len(speciess)):
		if genera[i]!=genus:
			continue
		
		xs.append(ess[i])
		ys.append(lowtemp_growths[i])
		subqvalues.append(qvalues[i])
		subspeciess.append(speciess[i])
	
	#line, = axis.plot(xs,ys,':',label=genus)
	#color = pylab.getp(line,'color')
	
	for i in range(0,len(xs)):
		
		if subqvalues[i]<0.05:
			symbol='o'
		else:
			symbol='s'
		
		if subspeciess[i]==species[1]:
			color='r'
		elif subspeciess[i]==species[2]:
			color='r'
		else:
			color='k'
		
		temp_axis.plot(xs[i],ys[i],symbol,markersize=3,color=color)
	
	#legend_axis.plot([-2,-1],[-2,-1],'.',color=color,markeredgewidth=0.0, label=genus)

temp_axis.plot([-1],[-1],'o',label='q<0.05',markersize=3,color='w',markeredgecolor='k')
temp_axis.plot([-1],[-1],'s',label='q>0.05',markersize=3,color='w',markeredgecolor='k')
temp_axis.legend(loc='upper right',frameon=False,numpoints=1,handletextpad=0.1)
host_axis.set_xlim([-0.05,1.05*max_host_d])
host_axis.set_ylim([-0.05,1.05*max_host_d])
host_hist_axis.set_ylim([-0.05,1.05*max_host_d])
host_hist_axis.set_xticklabels([])
host_hist_axis.set_yticklabels([])
host_hist_axis.set_xticks([])
host_host_hist_axis.set_xlim([-0.05,1.05*max_host_d])
host_host_hist_axis.set_xticklabels([])
host_host_hist_axis.set_yticklabels([])
host_host_hist_axis.set_yticks([])
for s in [1,2]:
	bacteria_axis[s].set_xlim([-0.05,1.05*max_host_d])
	bacteria_axis[s].set_ylim([-0.05,1.05*max_bacteria_d[s]])
	
	bacteria_axis[s].set_ylabel("%s phylogenetic distance" % get_pretty_name(species[s]))
	bacteria_axis[s].set_xlabel('Host phylogenetic distance')
	
	bacteria_hist_axis[s].set_ylim([-0.05,1.05*max_bacteria_d[s]])
	bacteria_hist_axis[s].set_xticklabels([])
	bacteria_hist_axis[s].set_yticklabels([])
	bacteria_hist_axis[s].set_xticks([])
	
	
	bacteria_host_hist_axis[s].set_xlim([-0.05,1.05*max_host_d])
	bacteria_host_hist_axis[s].set_xticklabels([])
	bacteria_host_hist_axis[s].set_yticklabels([])
	bacteria_host_hist_axis[s].set_yticks([])
	

bac_bac_axis.set_xlabel('%s phylogenetic distance' % get_pretty_name(species[2]))
bac_bac_axis.set_ylabel('%s phylogenetic distance' % get_pretty_name(species[1]))
bac_bac_axis.set_ylim([-0.05,1.05*max_bacteria_d[1]])
bac_bac_axis.set_xlim([-0.05,1.05*max_bacteria_d[2]])
bac_bac_hist_axis.set_ylim([-0.05,1.05*max_bacteria_d[1]])
bac_bac_hist_axis.set_xticklabels([])
bac_bac_hist_axis.set_yticklabels([])
bac_bac_hist_axis.set_xticks([])
bac_bac_bac_hist_axis.set_xlim([-0.05,1.05*max_bacteria_d[2]])
bac_bac_bac_hist_axis.set_xticklabels([])
bac_bac_bac_hist_axis.set_yticklabels([])	
bac_bac_bac_hist_axis.set_yticks([])	


#legend_axis.legend(loc='center',frameon=False,ncol=1,numpoints=1)
temp_axis.set_xlabel('PACo Effect Size')
temp_axis.set_ylabel('Relative growth at 27$^{\circ}$C')
temp_axis.set_ylim([-0.05,1.05])
temp_axis.set_xlim([-0.0075,0.0325])		
print("Saving fig!")
fig.savefig('%s_%s_distance_correlation.pdf' % (species[1],species[2]),bbox_inches='tight',transparent=True)
fig.savefig('%s_%s_distance_correlation.png' % (species[1],species[2]),bbox_inches='tight',dpi=300)
fig2.savefig('%s_%s_distance_cocorrelation.png' % (species[1],species[2]),bbox_inches='tight',dpi=300)
fig3.savefig('%s_%s_distance_within_country_cocorrelation.png' % (species[1],species[2]),bbox_inches='tight',dpi=300)

print("Done!")

