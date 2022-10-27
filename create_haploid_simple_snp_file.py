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
    sample_ids.append(record.id)
    snp_matrix.append(list(record.seq))
    
snp_matrix = numpy.array(snp_matrix)
snp_matrix = snp_matrix.T

new_snp_matrix = []

num_removed = 0
for i in range(0,snp_matrix.shape[0]):
	# Need to see if site is truly polymorphic. 
	gap_idxs = (snp_matrix[i,:]=='-')
	N_idxs = (snp_matrix[i,:]=='N')
	good_idxs = numpy.logical_not(numpy.logical_or(gap_idxs,N_idxs))
	
	n=good_idxs.sum() # sample size
	
	if n==0:
		num_removed+=1
		continue
		
	alleles = Counter(snp_matrix[i,good_idxs])
	if len(alleles)==1:
		# Monomorphic!
		num_removed+=1
		continue
	
	allele_counts = alleles.most_common()
	k = allele_counts[0][1]
	n = good_idxs.sum()
	MAF = (n-k)*1.0/n
	
	#if(alleles.most_common()[-1][1]<2):
		# Singleton! Don't really want that either?
		#continue
		
	if(MAF<0.1):
		# make it like human genome?
		continue

	snp_matrix[i,gap_idxs] = 'N'
		
	# A real snp
	#if(gap_idxs.sum()>0):
		#print("SNP with GAP!")
		#print(alleles)
		#print(snp_matrix[i,:])
		# Try to replace with N
		
	new_snp_matrix.append(snp_matrix[i,:])
	
snp_matrix = numpy.array(new_snp_matrix)

print("Num removed", num_removed)
sys.exit(0)

output_str = "\t".join(["Chrom", "Pos"]+sample_ids)
print(output_str)
for i in range(0,snp_matrix.shape[0]):
    if option=='1' and i>(snp_matrix.shape[0]/2):
        continue

    if option=='2' and i<(snp_matrix.shape[0]/2):
        continue

    genotypes = list(snp_matrix[i,:])
    
    output_str = "\t".join(["1",str(i+1)]+genotypes)
    print(output_str)

