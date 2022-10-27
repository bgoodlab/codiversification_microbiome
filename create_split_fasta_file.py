import sys
from Bio import SeqIO
import numpy
from collections import Counter

fasta_file = SeqIO.parse(sys.argv[1], "fasta")

if len(sys.argv) < 3:
    option = 'a'
else:
    option = sys.argv[2].strip()

sample_ids = []
snp_matrix = []
for record in fasta_file:
	print(">%s" % record.id)
	L = len(record.seq)
	if option=='1':
		print(str(record.seq)[:L//2])
	elif option=='2':
		print(str(record.seq)[L//2:])
	else:
		print(str(record.seq))
			
