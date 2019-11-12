#Part C - Adapting content from this course and your own independent learning to create new research (30 marks)

#For this part of the assignment I want to determine whether two genes from the same taxonomic group generate the same phylogenetic hypothesis. I will be using ribosomal RNA genes 16S and 18S to reconstruct phylogenetic trees for the phylum Daphnia. I wish to determine whether the same clusters or groupings of specific species can be observed in the trees generated for these genes.I will also try to determine which of these genes is better for inferring phylogenetic relationships for Daphnia. 

#######################################################################################################
#Load Packages
library(tidyverse)
library(rentrez)
library(muscle)
library(Biostrings)
library(DECIPHER)
library(dendextend)
library(ape)
library(phangorn)
library(phytools)
######################################################################################################

#Acquiring and cleaning up data for Daphnia 16S and  18S rRNA genes fron NCBI's GeneBank

#----------------------Processing and Filtering Data for Daphnia 16S rRNA gene------------------------

#Searching for Daphnia 16S rRNA specifically excluding full genomes
Daphnia16S_search <- entrez_search(db = "nuccore", term = "Daphnia[ORGN] AND 16S[TITL] NOT (genome[TITL])", retmax = 350)

#Attributes checks
Daphnia16S_search$ids
Daphnia16S_search$count

#Fetching data 
Daphnia16S_fetch <- entrez_fetch(db = "nuccore", id = Daphnia16S_search$ids, rettype = "fasta")

#Attribute checks
class(Daphnia16S_fetch)
#head(Daphnia16S_fetch)

#Downloading data in fasta format: 11:35 AM October 30,2019 
write(Daphnia16S_fetch, "Daphnia_fetch16.fasta", sep = "\n")

#------------------------------------Converting data to a data frame-----------------------------------------

#Extracting data from file 
stringset_D16 <- readDNAStringSet("Daphnia_fetch16.fasta")

#Converting data set into a data frame
dfD16 <- data.frame(Original = names(stringset_D16), Sequence = paste(stringset_D16))
#writing original data into file
write.csv(dfD16, "Daphnia16S_OGData.csv")

#Cleaning up names and separating them into new columns
dfD16$Species_Name <- word(dfD16$Original, 2L, 3L)
dfD16$Gene_Name <- word(dfD16$Original, 6L)
dfD16$Unique_identifier <- word(dfD16$Original, 1L)

#Rearranging columns
dfD16 <- dfD16[, c("Unique_identifier", "Species_Name","Gene_Name", "Sequence", "Original")]

#Number of total sequences
length(dfD16$Sequence)

#Number of unique species
length(unique(dfD16$Species_Name))

#Checking how many unqiue genes in data set
length(unique(dfD16$Gene_Name))

#Checking for mininum, max, median and mean sequence lengths
min(str_length(dfD16$Sequence))
mean(str_length(dfD16$Sequence))
median(str_length(dfD16$Sequence))
max(str_length(dfD16$Sequence))

#Checking distribution of sequence lengths
hist(str_length(dfD16$Sequence),xlab = 'Sequence Length', ylab = 'Frequency', main = 'Fig 1. Distribution of 16S rRNA Sequence Lengths Frequency')


#Filtering sequences longer than 600 bp
D16_filtered <- dfD16 %>%
  filter(!str_length(dfD16$Sequence) > 600)

#Checks to make sure >600 bp sequences are removed
min(str_length(D16_filtered$Sequence))
max(str_length(D16_filtered$Sequence))
median(str_length(D16_filtered$Sequence))
mean(str_length(D16_filtered$Sequence))

#Checking for the number of unique species after filtering
length(unique(D16_filtered$Species_Name))


#----------------------Processing and Filtering Data for Daphnia 18S rRNA gene------------------------

#Searching for Daphnia 16S rRNA specifically excluding full genomes
Daphnia18S_search <- entrez_search(db = "nuccore", term = "Daphnia[ORGN] AND 18S[TITL] NOT (genome[TITL])", retmax = 350)

#Attributes check
Daphnia18S_search$ids
Daphnia18S_search$count

#Fetching search data 
Daphnia18S_fetch <- entrez_fetch(db = "nuccore", id = Daphnia18S_search$ids, rettype = "fasta")

#Attribute checks
class(Daphnia18S_fetch)
#head(Daphnia16S_fetch)

#Downloading data in fasta format: 11:35 AM October 30,2019 
write(Daphnia18S_fetch, "Daphnia_fetch18.fasta", sep = "\n")

#------------------------------------Converting data to a data frame-----------------------------------------

#Extracting data from file 
stringset_D18 <- readDNAStringSet("Daphnia_fetch18.fasta")

#Converting data set into a data frame
dfD18 <- data.frame(Original = names(stringset_D18), Sequence = paste(stringset_D18))
#writing original data into file
write.csv(dfD18, "Daphnia16S_OGData.csv")

#Cleaning up names and separating them into new columns
dfD18$Species_Name <- word(dfD18$Original, 2L, 3L)
dfD18$Gene_Name <- word(dfD18$Original, 6L)
dfD18$Unique_identifier <- word(dfD18$Original, 1L)

#Rearranging columns
dfD18 <- dfD18[, c("Unique_identifier", "Species_Name","Gene_Name", "Sequence", "Original")]

#Number of total sequences
length(dfD18$Sequence)

#Number of unique species
length(unique(dfD18$Species_Name))

#Checking how many unqiue genes in data set
length(unique(dfD18$Gene_Name))

#Checking for mininum, max, median and mean sequence lengths
min(str_length(dfD18$Sequence))
mean(str_length(dfD18$Sequence))
median(str_length(dfD18$Sequence))
max(str_length(dfD18$Sequence))

#Checking distribution of sequence lengths
hist(str_length(dfD18$Sequence),xlab = 'Sequence Length', ylab = 'Frequency', main = 'Fig 2. Distribution of 18S rRNA Sequence Lengths Frequency')


#Filtering sequences longer than 600 bp
D18_filtered <- dfD18 %>%
  filter(!str_length(dfD18$Sequence) > 600)

D18_Clean <- D18_filtered %>%
  filter(!str_length(D18_filtered$Sequence) < 450)

#Checks to make sure >600 bp sequences are removed
min(str_length(D18_Clean$Sequence))
max(str_length(D18_Clean$Sequence))
median(str_length(D18_Clean$Sequence))
mean(str_length(D18_Clean$Sequence))


#Checking for the number of unique species after filtering
length(unique(D18_Clean$Species_Name))

#After cleaning up and filtering the data sets for both 16S and 18S rRNA gene, it is clear that more work has been done with 16S rRNA data. This may mean that the quality of data for the 16S will likely be higher than that of 18S. There are also a much higher number of unique species (higher recorded diversity) in the GeneBank for 16S. These factors may lead to a much better clustering and overall more reliable data for 16S. 

#####################################################################################################

#Multiple Sequence Alignments(MSA)

#----------------------------------- MSA for Daphnia 16S--------------------------------------------

#Further filtering out a sample that has way too many "N"s and possible polymorphisms - re-ran MSA 
D16_filtered <- D16_filtered[-320,]

#Reformat sequences into DNAStringSet
D16_filtered$Sequence <-DNAStringSet(D16_filtered$Sequence)

#assigning the species names to their nucleotides as labels throughout dowsntream analysis
names(D16_filtered$Sequence) <- D16_filtered$Species_Name

#Performing MSA using the function muscle from the muscle package. All options are left to default.
#I will be using the same parameters throughout the entire project for both 16S and 18S. I chose muscle for its reliability in performing MSA and for its speed. I did not want to reduce or increase the gap penalty as both rRNA genes are highly conserved. I want to allow insertions of gaps in account for point deletions but also do not want excessive insertions that may create false alignments. 
D16.alignment <- DNAStringSet(muscle::muscle(D16_filtered$Sequence), use.names = TRUE)


#Display multiple sequence alignment - showed one sequence to have poor quality (sequence #320)
#BrowseSeqs(D16.alignment)

#Calculating mean, max, min and observing distribution of gaps through a histogram
min(unlist(lapply(D16.alignment, str_count, "-")))
max(unlist(lapply(D16.alignment, str_count, "-")))
mean(unlist(lapply(D16.alignment, str_count, "-")))
hist(unlist(lapply(D16.alignment, str_count, "-")))
hist(str_length(dfD16$Sequence),xlab = 'Number of Gaps', ylab = 'Frequency', main = 'Fig 3. Distribution of Gaps in 16S sequence MSA')


#Converting data type into DNAbin (using as.DNAbin from ape package) for more downstream analysis
dnaBin.D16 <- as.DNAbin(D16.alignment)

#Creating a distance matrix using the "K80" model
#The model K80 was chosen due to its frequent use in literature and produces very similar distance matrices to other methods. Pairwise deletion is set to true therefore only missing values that are not found in all variables will be deleted. I don't want to lose too much data that may be relevant, especially since the data sets I am using are not very large. Again, the same parameters for calculating distance for 16S MSA will be applied to 18S as well. 
DM.D16 <- dist.dna(dnaBin.D16, model = "K80", as.matrix = TRUE, pairwise.deletion = TRUE)

#Clustering using the UPGMA method with a 2% divergence threshold 
#I will be using UPGMA(or average-linkage) as my method of choice for clustering all the DNA distances for both rRNA data throughout. UPGMA was chosen because it improves upon the limitations of single-linkage (connecting distant relatives closely) and complete-linkage (compact clustering).
clusters.D16 <- IdClusters(DM.D16,
                           method = "UPGMA",
                           cutoff= 0.02,      #2% sequence difference as species boundary
                           showPlot = TRUE,   #returns a dendogram
                           type = "clusters", #returns a clusters output
                           verbose = TRUE)
title("Fig 4. Daphnia Phylogeny Based on 16S rRNA")

#----------------------------------- MSA for Daphnia 18S--------------------------------------------


#Reformat sequences into DNAStringSet
D18_Clean$Sequence <-DNAStringSet(D18_Clean$Sequence)

#Assigning the species names to their nucleotides as labels throughout dowsntream analysis
names(D18_Clean$Sequence) <- D18_Clean$Species_Name

#Performing MSA using the function muscle from the muscle package. All options are left to default.
D18.alignment <- DNAStringSet(muscle::muscle(D18_Clean$Sequence), use.names = TRUE)


#Display multiple sequence alignment in a new internet browser window
BrowseSeqs(D18.alignment)

#Calculating mean, max, min and observing distribution of gaps through a histogram
min(unlist(lapply(D18.alignment, str_count, "-")))
max(unlist(lapply(D18.alignment, str_count, "-")))
mean(unlist(lapply(D18.alignment, str_count, "-")))
hist(unlist(lapply(D18.alignment, str_count, "-")))
hist(str_length(dfD16$Sequence),xlab = 'Number of Gaps', ylab = 'Frequency', main = 'Fig 5. Distribution of Gaps in 16S sequence MSA')


#Converting data type into DNAbin for more downstream analysis
dnaBin.D18 <- as.DNAbin(D18.alignment)

#Creating a distance matrix using the "K80" model
DM.D18 <- dist.dna(dnaBin.D18, model = "K80", as.matrix = TRUE, pairwise.deletion = TRUE)

#Clustering using the UPGMA method with a 2% divergence threshold 
clusters.D18 <- IdClusters(DM.D18,
                           method = "UPGMA",
                           cutoff= 0.02,      #2% sequence difference as species boundary
                           showPlot = TRUE,   #Generates a dendogram
                           type = "clusters", #returns a clusters output
                           verbose = TRUE)
title("Fig 6. Daphnia Phylogeny Based on 18S rRNA")

#The MSA for 18S produced much higher number of gaps compared to the MSA for 16S. This means a lot more gaps are inserted in order to create the best alignments. This may mean that the 18S gene sequences for the species acquired for this project have a lot of differences in their gene sequences caused by deletions. 

#Graph is made up of distinct clusters that have very high dissimilarities and has a very large distance scale. The DNA distances are very low within clusters but are extremely large in between clusters. Sequences of records (such as that of data 1) that seems to be outliers we're searched in BLASTn to check if they were mislabled. The searches clarified that these sequences are actually Daphnia 18S rRNA sequences and matched with sequences from the same species. 

##############################################################################################################

#------------Extracting one record per species from both Daphnia 16S and 18S rRNA datasets--------------------

#The dendograms generated by including all the sequence are too crowded and are not easy to visually analyze.
#Randomly selecting one species might bring a much better visualization and give a more conclusive suggestion on the phylogenetic relationships of Daphnia.

#-------------------Distance matrix and re-run of MSA for 16S using representative species------------

#Random selection of one representative sequence per species
D16_SeqPer_Species <- D16_filtered %>%
  group_by(Species_Name) %>%
  sample_n(1)


#Checking the number of unique species
length(unique(D16_SeqPer_Species$Species_Name))


#Converting sequences into a DNAStringSet data type
D16_SeqPer_Species$Sequence <- DNAStringSet(D16_SeqPer_Species$Sequence)

#Setting species name as identifiers for each sequence
names(D16_SeqPer_Species$Sequence) <- D16_SeqPer_Species$Species_Name


#Performing multiple sequence alignment of each representative sequence from each species
D16_SeqPer_Species.alignment <- DNAStringSet(muscle::muscle(D16_SeqPer_Species$Sequence), use.names = TRUE)


#Converting alignment data into a DNAbin
D16.Per.species.DNAbin <- as.DNAbin(D16_SeqPer_Species.alignment)

#Creating a distance matrix using the "K80" model
DM.D16.Per.species <- dist.dna(D16.Per.species.DNAbin, model = "K80", as.matrix = TRUE, pairwise.deletion = TRUE)


#-------------------Distance matrix and re-run of MSA for 16S using representative species-------------

#Random selection of one representative sequence per species
D18_SeqPer_Species <- D18_Clean %>%
  group_by(Species_Name) %>%
  sample_n(1)


#Checking the number of unique species
length(unique(D18_SeqPer_Species$Species_Name))


#Converting sequences into a DNAStringSet data type
D18_SeqPer_Species$Sequence <- DNAStringSet(D18_SeqPer_Species$Sequence)

#Setting species name as identifiers for each sequence
names(D18_SeqPer_Species$Sequence) <- D18_SeqPer_Species$Species_Name


#Performing multiple sequence alignment of each representative sequence from each species
D18_SeqPer_Species.alignment <- DNAStringSet(muscle::muscle(D18_SeqPer_Species$Sequence), use.names = TRUE)


#Converting alignment data into a DNAbin
D18.Per.species.DNAbin <- as.DNAbin(D18_SeqPer_Species.alignment)

#Creating a distance matrix using the "K80" model
DM.D18.Per.species <- dist.dna(D18.Per.species.DNAbin, model = "K80", as.matrix = TRUE, pairwise.deletion = TRUE)

###############################################################################################################


#-------------------------Visualizing the 18S and 16S phylogenetic trees individually--------------------------

#Creating a dendogram for 16S using the function hclust to create a cluster and setting dendogram parameters using denextend's set(den..) function.
tree_D16S <- D16.Per.species.DNAbin %>% 
  dist %>% 
  hclust(method="average")%>% #clustering based on average-linkage methodology 
  as.dendrogram %>% #converting cluster to a dendrogram
  set("branches_k_color", k=11) %>% #setting color based on clusters for branches
  set("labels_cex", c(.6,.6)) %>% #resizing label size
  set("nodes_pch", 19) %>% #shapes to mark nodes
  set("nodes_col", c("black"))%>% #color of nodes
  set("hang_leaves") #creates a hanging dendogram 
plot(tree_D16S,main = "Fig 7. UPGMA Tree: Daphnia 16S rRNA" )

#The phylogenetic tree generated for Daphnia 16S (Fig 7) is easier to read than the dendrogram which include all the sequences. This tree displays 7 unique groupings of species that are closely related to eachother. 

#Creating a dendogram for 18S using the function hclust to create a cluster and setting dendrogram parameters using denextend's set(den..) function.
tree_D18S <- D18.Per.species.DNAbin %>% 
  dist %>% 
  hclust(method="average")%>% #clustering based on average-linkage methodology 
  as.dendrogram %>% #converting cluster to a dendrogram
  set("branches_k_color", k=4) %>% #setting color based on clusters for branches
  set("labels_cex", c(.6,.6)) %>% #resizing label size
  set("nodes_pch", 19) %>% #shapes to mark nodes
  set("nodes_col", c("black"))%>% #color of nodes
  set("hang_leaves") #creates a hanging dendogram 
plot(tree_D18S, main = "Fig 8. UPGMA Tree: Daphnia 18S rRNA")

#The phylogenetic tree generated for Daphnia 18S (Fig 8) looks much cleaner than the previous dendogram. This tree very clearly displays 4 unique groupings of species that are closely related to eachother. 

#The distance scales on both of these trees are omitted as I am only interested on whether the same groupings of species can be observed in both graphs. 

##############################################################################################################

#---------------------Comparing the dendograms generated by using 18S and 16S rRNA genes----------------------

#I want to compare the two trees in the same image, and the function tanglegram from the package dendextend allows mirror comparisons of trees. I am aware that there is a extremely large discrepancy in the number of species for 18S ans 16S. Luckily, tanglegram prunes the the branches that the two trees do not share and only generate a comparison of shared variables.

#Creating a dendogram using 16S MSA derived distance matrix using hclust (average-linkage) from the base package
dd_D16S <- D16.Per.species.DNAbin %>% 
  dist %>% 
  hclust(method="average")%>% #clustering based on average-linkage methodology 
  as.dendrogram %>% #converting cluster to a dendrogram
  set("branches_k_color", k=11)  #setting color based on clusters for branches

#Creating a dendrogram using 18S MSA derived distance matrix using hclust (average-linkage) from the base package
dd_D18S <- D18.Per.species.DNAbin %>% 
  dist %>% 
  hclust(method="average")%>% #clustering based on average-linkage methodology 
  as.dendrogram %>% #converting cluster to a dendrogram
  set("branches_k_color", k=4) #setting color based on clusters for branches

#Compares the dendograms constructed from 18S rRNA and 16S rRNA using tanglegram from then dendextend package
tanglegram(dd_D16S,dd_D18S,main ="Fig 9. Phylogeny Reconstruction by rRNA", main_left="Daphnia 18S rRNA", main_right = "Daphnia 16S rRNA", common_subtrees_color_branches = TRUE, columns_width = c(5, 2, 5),margin_inner = 11,cex_main = 1, cex_main_left = 1, cex_main_right = 1, sort = FALSE, rank_branches = TRUE)

############################################################################################################

#----------------------------------------------Conclusions---------------------------------------------------

#In conclusion, the 16S rRNA gene and 18S rRNA gene do not produce the same phylogenetic hypothesis. In Fig 9, it can be observed that some species are placed in different clusters. For example, Daphnia curvirostris and Daphnia galeata were put in the same node by the of 18S rRNA tree but they were in the different nodes in the 16S tree. In Fig 4 and 6, we can we see that the distance scales in both clusters are drastically different. The distance scale for the 18S dendrogram was 4 times that of 16S. This means the clusters are so drastically different from each other and may suggest that making a tree based on 18S is not viable. The difference in the dendogram structure for both 16S and 18S and may have been due to the difference of available data. The clustering and DNA distance methods used for the study may also have not been the best fit for the data sets. A lot of improvements can be implemented to improve the results of this mini-project. The first would to perform a cophenetic correlation test on all the dendograms generated to determine their validity. Trying all the other available methods and performing cophoenetic correlation test on them will give an insight into which methods might have been the best choice. Furthermore, the data used for this mini-project only came from one source. So pulling other Daphnia 18S and 16S data from other sources and including them in the study will likely to improve the results of this study. 
