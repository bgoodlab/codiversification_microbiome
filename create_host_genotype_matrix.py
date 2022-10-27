import sys
from Bio import SeqIO
import numpy

id_id_map = {}
file = open("leylabdata/ID_table_n839.txt","r")
file.readline() # header
for line in file:
	items = line.split()
	fasta_id = items[1].strip()
	tree_id = items[0].strip()
	
	id_id_map[fasta_id]=tree_id
file.close()

genotype_conversion_map = {'A':'AA','T':'TT','G':'GG','C':'CC','N':'NN','R':'AG','Y':'CT','S':'GC','W':'AT','K':'GT','M':'AC'}
sample_ids = []
snp_matrix = []
fasta_file = SeqIO.parse("leylabdata/HostTree_alignment.fasta", "fasta")

regular_alleles = set(['A','C','T','G'])

if len(sys.argv) < 2:
    option = 'a'
else:
    option = sys.argv[1].strip()

for record in fasta_file:
    sample_ids.append(id_id_map[record.id])
    snp_matrix.append(list(record.seq))
    
snp_matrix = numpy.array(snp_matrix)
snp_matrix = snp_matrix.T

genotypess = []
for i in range(0,snp_matrix.shape[0]):
    if option=='1' and i>(snp_matrix.shape[0]/2):
        continue

    if option=='2' and i<(snp_matrix.shape[0]/2):
        continue

    ancestral_allele = ""
    derived_allele = ""
    allele_key = {}

    genotypes = []
    for g in snp_matrix[i,:]:
        if g in regular_alleles:
            if ancestral_allele=="":
                ancestral_allele=g
                allele_key[ancestral_allele]=0

            if g not in allele_key:
                if derived_allele=="":
                    derived_allele=g
                    allele_key[derived_allele]=1
                else:
                    print("ERROR!", g, allele_key, derived_allele,ancestral_allele)
                    sys.exit(1)
            genotypes.append(allele_key[g])
        elif g=='N':
            genotypes.append(-1)
        else:
            genotypes.append(0.5)

    genotypess.append(genotypes)
    
print(", ".join(sample_ids))
for i in range(0,len(genotypess)):
    print(", ".join([str(g) for g in genotypess[i]]))
