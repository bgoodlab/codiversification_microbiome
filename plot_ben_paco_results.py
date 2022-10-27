import pylab
import sys
import numpy
import matplotlib.colors as colors
import matplotlib.cm as cmx
from math import log10,ceil
import matplotlib as mpl
import matplotlib.gridspec as gridspec
from numpy.random import randint

mpl.rcParams['font.size'] = 7
mpl.rcParams['lines.linewidth'] = 0.5
mpl.rcParams['legend.frameon']	= False
mpl.rcParams['legend.fontsize']	 = 'small'

fig = pylab.figure(figsize=(5, 2))

# Set up grids to hold figure panels
outer_grid = gridspec.GridSpec(2, 1, height_ratios=[1,1], hspace=0.05)

global_axis = pylab.Subplot(fig, outer_grid[0])
fig.add_subplot(global_axis)

country_axis = pylab.Subplot(fig, outer_grid[1])
fig.add_subplot(country_axis)

def parse_ben_paco_results(filename):
	file = open(filename,"r")
	file.readline()# blank line at top
	
	speciess = []
	pvalues = []
	effect_sizes = []
	host_effect_sizes = []
	line = file.readline()
	while line!="":
		species = line.split(".")[1][3:]
		
		line = file.readline() # next one is pvalue
		if(not line.startswith('[1] "p-value')):
			print("Bad species: %s" % species)
			continue
		pvalue = float(line.strip().split("=")[1][:-1])
		
		
		line = file.readline() # next one is m2obs
		m2obs = float(line.strip().split("=")[1][:-1])
		
		line = file.readline() # next one is m2null
		m2null = float(line.strip().split("=")[1][:-1])
		
		line = file.readline() # next one is standard deviation
		pass
		
		line = file.readline() # next one is effect size
		effect_size = float(line.strip().split("=")[1][:-1])
		
		line = file.readline() # next one is host-host pvalue
		file.readline() # next one is host-host m2obs
		file.readline() # next one is host-host m2null
		file.readline() # next one is host-host m2std
		
		line = file.readline() # next one is host-host effect size
		host_effect_size = float(line.strip().split("=")[1][:-1])
		
		
		speciess.append(species)
		pvalues.append(pvalue)
		effect_sizes.append(effect_size)
		host_effect_sizes.append(host_effect_size)
		
		line = file.readline() # next one is next species
		
	return speciess,pvalues,effect_sizes, host_effect_sizes
	
def calculate_qvalues(pvalues):

	# Uses the formula
	#
	# Qi = min_{Q>Pi} { Q*(sum_j 1)/(sum_j theta(Q-Pj)) } 
	#
	# Qi = max_{P_(k)>Pi} 

	# calculate q-values
	qvalues = []
	
	sorted_pvalues = numpy.array(sorted(pvalues))

	Ntot = len(sorted_pvalues)
	Nless = numpy.array([(sorted_pvalues<=p).sum() for p in sorted_pvalues])

	for p in pvalues:

		min_q = 1e06
			
		for j in reversed(range(0,Ntot)):
			 
			if sorted_pvalues[j]<p:
				break
			else:
				new_q = Ntot*sorted_pvalues[j]*1.0/Nless[j]
				if new_q<min_q:
					min_q = new_q
					
		qvalues.append(min_q)
	
	qvalues = numpy.array(qvalues)
	
	return qvalues
	
speciess,pvalues,ess,host_ess = parse_ben_paco_results("paco_output_all_1000_nocluster.txt")
dummy1,country_pvalues,country_ess,country_host_ess = parse_ben_paco_results("paco_output_country_1000.txt")


qvalues = calculate_qvalues(pvalues)
country_qvalues = calculate_qvalues(country_pvalues)


ess, speciess,qvalues,host_ess, country_qvalues, country_ess, country_host_ess = zip(*sorted(zip(ess, speciess,qvalues,host_ess, country_qvalues, country_ess, country_host_ess),reverse=True))
speciess = list(speciess)

# Shouldn't matter, but necessary for plotting log scale
ess = numpy.fabs(ess)
host_ess = numpy.fabs(host_ess)
country_ess = numpy.fabs(country_ess)
country_host_ess = numpy.fabs(country_host_ess)

xs = range(0,len(speciess))
	
for desired_axis, desired_qvalues, desired_ess, desired_host_ess in zip([global_axis, country_axis], [qvalues, country_qvalues],[ess, country_ess],[host_ess,country_host_ess]):
	for i in range(0,len(speciess)):
		if desired_qvalues[i]<0.05:
			#print(speciess[i])
			#speciess[i] = ("**"+speciess[i])
			color='r'
		else:
			color='k'
		#pylab.plot([xs[i],xs[i]],[0,ess[i]],'-',color=color)
		desired_axis.plot([xs[i]],[desired_ess[i]],'o',color=color,markersize=2)
		desired_axis.plot([xs[i]],[desired_host_ess[i]],'v',color=color,markersize=1.5)
		desired_axis.plot([xs[i], xs[i]], [desired_ess[i],desired_host_ess[i]], ':', color=color, linewidth=0.5)
	
global_axis.plot([xs[-1]+100],[1],'o',color='k',markersize=2,label='Host vs Bacteria')
global_axis.plot([xs[-1]+100],[1],'v',color='k',markersize=2,label='Host vs Host')
global_axis.plot([xs[-1]+100],[1],'s',color='r',markersize=2,label='q<0.05')
global_axis.legend(loc='lower left',frameon=False,ncol=3,numpoints=1,handletextpad=0.2)
global_axis.set_xticks(xs)
#global_axis.set_xticklabels(speciess,rotation=90,fontsize=5)
global_axis.set_ylabel('Global ES')
#global_axis.set_ylim([2e-04,1e-01])
#global_axis.set_yticks([1e-03,1e-02,1e-01,1])
global_axis.set_xlim([-1,xs[-1]+1])
global_axis.spines['top'].set_visible(False)
global_axis.spines['right'].set_visible(False)

country_axis.set_xticks(xs)
country_axis.set_xticklabels(speciess,rotation=90,fontsize=5)
country_axis.set_ylabel('W/in country ES')
country_axis.set_ylim([2e-04,1])
country_axis.set_yticks([1e-03,1e-02,1e-01,1])

country_axis.set_xlim([-1,xs[-1]+1])
country_axis.spines['top'].set_visible(False)
country_axis.spines['right'].set_visible(False)


#pylab.gca().set_yticks([1e-03,1e-02,1e-01,1])
pylab.savefig('ben_paco_results.pdf',bbox_inches='tight')





