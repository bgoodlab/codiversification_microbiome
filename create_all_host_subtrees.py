import sys
import os

speciess = []

file = open("paco_results.csv","r")
for line in file:
    items = line.split(",")
    species = items[0].strip()

    os.system('python create_host_subtree.py %s' % species)
