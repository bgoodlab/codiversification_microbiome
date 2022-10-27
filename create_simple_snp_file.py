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

if len(sys.argv) < 2:
    option = 'a'
else:
    option = sys.argv[1].strip()

for record in fasta_file:
    sample_ids.append(id_id_map[record.id])
    snp_matrix.append(list(record.seq))
    
snp_matrix = numpy.array(snp_matrix)
snp_matrix = snp_matrix.T


output_str = "\t".join(["rs#", "alleles", "chrom", "pos", "strand", "assembly#", "center", "protLSID", "assayLSID", "panelLSID", "QCcode"]+sample_ids)
print(output_str)
for i in range(0,snp_matrix.shape[0]):
    if option=='1' and i>(snp_matrix.shape[0]/2):
        continue

    if option=='2' and i<(snp_matrix.shape[0]/2):
        continue

    genotypes = list(genotype_conversion_map[item] for item in snp_matrix[i,:])
    alleles = set("".join(genotypes))
    if 'N' in alleles:
        alleles.remove('N')
    allele_str = "/".join(list(alleles))
    output_str = "\t".join(["rs%d" % (i+1),allele_str,"1",str(i+1),"."]+(["NA"]*6)+genotypes)
    print(output_str)

