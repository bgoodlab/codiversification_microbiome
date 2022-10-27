import matplotlib  
import pylab
import sys
import numpy
from math import log10, fabs, log

import matplotlib.colors as colors
import matplotlib.cm as cmx
from math import log10,ceil
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from numpy.random import randint





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
	species = items[0].strip()
	es = float(items[4])
	phylum = items[3].strip()
	lowtemp_growth = float(items[6])
	qvalue = float(items[7])
	
	#genus = (species.split("_")[0]).split()[0]
	genus = items[8].strip() # actually family
	speciess.append(species)
	ess.append(es)
	lowtemp_growths.append(lowtemp_growth)
	phyla.append(phylum)
	genera.append(genus)
	qvalues.append(qvalue)
	
	
allowed_genera = set(genera)

NUM_COLORS = len(allowed_genera)

cm = pylab.get_cmap('jet')
#mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=[cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)]) 
mpl.rcParams['font.size'] = 7
mpl.rcParams['lines.linewidth'] = 0.5
mpl.rcParams['legend.frameon']  = False
mpl.rcParams['legend.fontsize']  = 'small'

fig = plt.figure(figsize=(5, 2.5))

# Set up grids to hold figure panels
outer_grid = gridspec.GridSpec(1, 2, width_ratios=[1,0.5], wspace=0.1)

axis = plt.Subplot(fig, outer_grid[0])
fig.add_subplot(axis)

legend_axis = plt.Subplot(fig, outer_grid[1])
fig.add_subplot(legend_axis)

legend_axis.set_ylim([0,1])
legend_axis.set_xlim([0,1])

legend_axis.spines['top'].set_visible(False)
legend_axis.spines['right'].set_visible(False)
legend_axis.spines['left'].set_visible(False)
legend_axis.spines['bottom'].set_visible(False)

legend_axis.set_xticks([])
legend_axis.set_yticks([])  

for genus in sorted(allowed_genera):
	xs = []
	ys = []
	subqvalues = []
	for i in range(0,len(speciess)):
		if genera[i]!=genus:
			continue
		
		xs.append(ess[i])
		ys.append(lowtemp_growths[i])
		subqvalues.append(qvalues[i])
	
	line, = axis.plot(xs,ys,':',label=genus)
	color = pylab.getp(line,'color')
	
	for i in range(0,len(xs)):
		
		if subqvalues[i]<0.05:
			symbol='o'
		else:
			symbol='s'
			
		axis.plot(xs[i], ys[i],symbol,label=genus,markeredgecolor='k',markeredgewidth=0.25,markersize=3,color=color)
	
	legend_axis.plot([-2,-1],[-2,-1],'.',color=color,markeredgewidth=0.0, label=genus)
	

legend_axis.legend(loc='center',frameon=False,ncol=1,numpoints=1)
axis.set_xlabel('PACo Effect Size')
axis.set_ylabel('Growth at 27C')
fig.savefig('temperature_correlation.pdf',bbox_inches='tight')