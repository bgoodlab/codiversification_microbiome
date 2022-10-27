####
#
# This code is lightly modified from the PACo script provided by Suzuki et al (Science 2022)
#
####
library(ape)
library(phylotools)
library(phytools)
library(ggtree)
library(dynamicTreeCut)
library(stringr)
library(cluster)
# stuff I needed to add:
library(reshape2)
library(vegan)

##2. Paco function
#This function takes 5 inputs and calculate the Paco test statistics and p-values using N permutations
#1. Bac_name: name of bacteria
#2. Bac_tree: bacterial tree (best tree or consensus tree)
#3. Host_tree: host tree (best tree or consensus tree)
#4. N_perm: the number of permutations for a given cospeciation function
#5. output_dir: output directory of the resulting file 

PACO.cutree.function = function(Bac_name, Bac_tree, Host_tree, N_perm, host_country_map, use_clusters){

#RemoveDups to select unique bacterial strain to create cuttree 
RemoveDups <- function(df, column) {
  inds = sample(1:nrow(df))  
  df   = df[inds, ]

  dups = duplicated(df[, column])
  df   = df[!dups, ]
  inds = inds[!dups]

  df[sort(inds, index=T)$ix, ]
}


print(use_clusters)
##Main function  
##1. Create Bac tree using cutree
#Get tip labels
Bac_tip_label = Bac_tree$tip.label

#Convert to distance matrix
Bac_dist = cophenetic(Bac_tree) 
#print("Done with cophenetic")
#Convert to hclust
Bac_hclust = diana(Bac_dist)
#print("Done with diana")
Bac_hclust_1 = as.hclust(Bac_hclust)
#print("Done with hclust")
 
#Cutree
Bac_cutree = cutreeDynamic(Bac_hclust_1, cutHeight = 0.99, method = "tree", minClusterSize = 2) # vector of cluster labels
#print("Bac_cutree structure")
#print(str(Bac_cutree))
Bac_cutree2 = cbind(Bac_tip_label,Bac_cutree) # create n x 2 matrix of tip label vs cluster label
#print("Bac_cutree2 structure")
#print(Bac_cutree2[1,])

Bac_cutree3 = as.data.frame(Bac_cutree2) # Turn into data frame
Bac_cutree3$Bac_cutree = sub("^", "X", Bac_cutree3$Bac_cutree) # Ad X at beginning of cluster labels (don't really know why)
#print("Bac_cutree3 structure")
#print(str(Bac_cutree3))

#Filter cutree IDs
Bac_cutree3_unique = RemoveDups(Bac_cutree3, "Bac_cutree") # remove only keep one strain per cluster
Bac_cutree3_unique_ID = Bac_cutree3_unique$Bac_tip_label # tip labels for those unique ones

#Filter bac tree
BacTree_cutree = keep.tip(Bac_tree, Bac_cutree3_unique_ID) # only keep unique IDs

#rename bac tree
BacTree_cutree_renamed = sub.taxa.label(BacTree_cutree, Bac_cutree3_unique) # added unique relative to the original version. 
Bac_dist_cutree = cophenetic(BacTree_cutree_renamed) 

##2. Edit host tree
#Host tree with IDs that exist in Bac tree
Host_tree_filter = keep.tip(Host_tree, Bac_tip_label)
Host_tree_filter_m = cophenetic(Host_tree_filter)

##3. Create HP file
Bac_cutree3_filtered = dcast(Bac_cutree3, Bac_tip_label~Bac_cutree, length, value.var = "Bac_cutree")
Bac_cutree3_filtered2 = data.frame(Bac_cutree3_filtered[,-1], row.names = Bac_cutree3_filtered[,1])

host_dist = Host_tree_filter_m #host tree filtered by bac tree IDs
host_names = colnames(host_dist)
nhosts = length(host_names)
#print("Done with cuttree stuff")		
#4. Rename variables to input to PACO function
if(use_clusters){
	Bac_dist = Bac_dist_cutree #bac tree that are representing one strain per cutree groups
	bac_names = colnames(Bac_dist)
	
	HP_file = as.matrix(Bac_cutree3_filtered2)
	HP = HP_file
} else{
	Bac_dist = Bac_dist 
	bac_names = colnames(Bac_dist)
	
	HP_file = matrix(0, nhosts,nhosts)

	for (i in 1:nhosts){
		for (j in 1:nhosts){
			if(host_names[i]==bac_names[j]){
				HP_file[i,j]=1
			}
		}
	}
	HP = HP_file
}

summary(warnings())

#5. PACO function
HP.bin <- which(HP > 0, arr.in=TRUE)
H.PCo <- pcoa(host_dist, correction="cailliez")$vectors #Performs PCo of Host distances 
P.PCo <- pcoa(Bac_dist, correction="cailliez")$vectors #Performs PCo of Parasite distances
H.PCo <- H.PCo[HP.bin[,1],] #adjust Host PCo vectors 
P.PCo <- P.PCo[HP.bin[,2],]  ##adjust Parasite PCo vectors
host_names = host_names[HP.bin[,1]]

# Now make country host map
country_host_map = list()
for(i in c(1:nhosts)){
	#print("Host name")
	#print(host_names[i])
	country = host_country_map[[host_names[i]]]
	#print(country)
	if(is.null(country_host_map[[country]])){
		#print("adding!")
		country_host_map[[country]] = c(i)
	}
	else{
	    country_host_map[[country]] = append(country_host_map[[country]], i) 
	}
	#print("Adding element")
	
}
#print(country_host_map)


HP.proc <- procrustes(H.PCo, P.PCo) #Procrustes Ordination 

m2.obs <- HP.proc$ss #observed sum of squares
N.perm = N_perm #set number of permutations for testing
P.value = 1
seed <-.Random.seed[trunc(runif(1,1,626))]
set.seed(seed)
    #set.seed(5) ### use this option to obtain reproducible randomizations

m2_null = c()

for (n in c(1:N.perm))
	{ 
	# do country specific permutation
	for(country in names(country_host_map)){
		original_idxs = country_host_map[[country]]
		if(length(original_idxs)>1){
			permuted_idxs = sample(country_host_map[[country]])
			H.PCo[original_idxs,] = H.PCo[permuted_idxs,]
		}
	}
	
	m2.perm <- procrustes(H.PCo, P.PCo)$ss #randomized sum of squares
	m2_null = append(m2_null,m2.perm)	
	if (m2.perm <= m2.obs)
		{P.value = P.value + 1} 
		
}
P.value <- P.value/(N.perm+1)

#Calculate m2_null mean and sd
mean_null = mean(m2_null)
sd_null = sd(m2_null)
print(paste0("p-value = ",P.value))
print(paste0("m2_obs = ",m2.obs))
print(paste0("Mean m2_null = ",mean_null))
print(paste0("SD m2_null = ",sd_null))
print(paste0("ES = ",(mean_null-m2.obs)/mean_null))
summary(warnings())
}

# Example
#Name of the bacteria to test
args = commandArgs(trailingOnly=TRUE)
Bac_filename = args[1]
Host_filename = args[2]
Host_1_filename = args[3]
Host_2_filename = args[4]

num_permutations = strtoi(args[5])
permutation_type = args[6] # allowed types = "--permutation=all", "--permutation=country", "--permutation=continent"
cluster_type = args[7] # allowed types = "--cluster=True", "--cluster=False"
if(cluster_type=='--cluster=True'){
	use_clusters = TRUE
} else {
	use_clusters = FALSE
}

print(Bac_filename)

#Best maximum likelihood tree out of 100 trees (output of StrainPhlAn)
Bac_tree = read.tree(Bac_filename) 
Bac_tip_label = Bac_tree$tip.label

#Best maximum likelihood tree out of 100 trees (output of SNPhylo)
Host_tree = read.tree(Host_filename)
Host_tree = keep.tip(Host_tree, Bac_tip_label)

Host_tree_1 = read.tree(Host_1_filename)
Host_tree_1 = keep.tip(Host_tree_1, Bac_tip_label)

Host_tree_2 = read.tree(Host_2_filename)
Host_tree_2 = keep.tip(Host_tree_2, Bac_tip_label)
# Create host country map

host_country_file = read.table(file = "leylabdata/ID_table_n839.txt", sep ="\t", header = TRUE)
num_hosts = dim(host_country_file)[1]
host_country_map = list()
for(i in c(1:num_hosts)){
	if(permutation_type=='--permutation=country'){
		host_country_map[[host_country_file[i,1]]] = host_country_file[i,6] # country
	} else if(permutation_type=='--permutation=continent'){
		host_country_map[[host_country_file[i,1]]] = host_country_file[i,7] # continent
	} else{
		host_country_map[[host_country_file[i,1]]] = 'all' # all
	}
}

PACO.cutree.function(Bac_name, Bac_tree, Host_tree, num_permutations, host_country_map, use_clusters)

PACO.cutree.function('host', Host_tree_1, Host_tree_2, 1000, host_country_map, FALSE)

