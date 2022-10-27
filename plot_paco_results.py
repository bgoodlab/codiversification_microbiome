import pylab
import numpy
import sys

speciess = []
ess = []
qvalues = []

file = open("paco_results.csv","r")
for line in file:
    items = line.split(",")
    species = items[0].strip()
    es = float(items[1])
    qvalue = float(items[2])

    speciess.append(species)
    ess.append(es)
    qvalues.append(qvalue)


ess, speciess,qvalues = zip(*sorted(zip(ess, speciess,qvalues),reverse=True))
speciess = list(speciess)

pylab.figure(figsize=(7,1))
xs = range(0,len(speciess))
for i in range(0,len(speciess)):
    if qvalues[i]<0.05:
        #print(speciess[i])
        #speciess[i] = ("**"+speciess[i])
        color='r'
    else:
        color='k'
    pylab.plot([xs[i],xs[i]],[0,ess[i]],'-',color=color)
    pylab.plot([xs[i]],[ess[i]],'.',color=color)

pylab.gca().set_xticks(xs)
pylab.gca().set_xticklabels(speciess,rotation=90,fontsize=6)
pylab.gca().set_ylabel('Effect size')
pylab.gca().set_ylim([0,0.1])
pylab.plot(xs,numpy.zeros_like(xs),'k:',linewidth=0.5)
pylab.gca().spines['top'].set_visible(False)
pylab.gca().spines['right'].set_visible(False)
pylab.savefig('paco_results.pdf',bbox_inches='tight')





