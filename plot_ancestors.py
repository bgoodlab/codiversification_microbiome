import numpy
from numpy.random import randint
import pylab
import sys
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

N = 10000

tmax = int(numpy.log2(N)*4)

ts = numpy.arange(0,tmax)

ns = [1]


for t in ts[1:]:
	
	parent1s = randint(0,N,size=ns[-1])
	parent2s = randint(0,N-1,size=ns[-1])
	parent_2s = parent2s+(parent2s>=parent1s)
	
	all_parents = numpy.hstack([parent1s,parent2s])
	
	n = len(set(all_parents))	
	#print n
	ns.append(n)
	
tgen = 25

tstar = numpy.log2(N)*tgen

print(tstar)

ts = ts*tgen
pylab.figure(1,figsize=(2,1))
#pylab.semilogy(ts*1.25,numpy.ones_like(ns)*N,'k:')
pylab.semilogy(ts,ns,'-',label='Autosomes')
pylab.semilogy(ts,numpy.ones_like(ns),'-',label='Mitochondria')
pylab.semilogy([tstar,tstar],[3e-01,3*N],'k:')
pylab.gca().set_ylim([3e-01,3*N])
pylab.gca().set_yticks([1,10,100,1000,10000])
pylab.gca().set_ylabel('# ancestors')
pylab.gca().set_xlabel('Years before present')
pylab.legend(frameon=False,loc='center')

labels = [item.get_text() for item in pylab.gca().get_yticklabels()]
labels[-1] = 'N'

pylab.gca().set_yticklabels(labels)

#pylab.gca().set_yticklabels(['1','','','','N'])
pylab.savefig('ancestors.pdf',bbox_inches='tight')
	