#Edited by Maverick Aradanas for Assignment 3
#All added comments will have a prefix of "MA"

####LOAD NECESSARY LIBRARIES####
library(rentrez)
library(seqinr)
library(Biostrings)
library(stringr)
library(tidyverse)
library(RSQLite)
library(DECIPHER)
library(muscle)
library(msa)
library(ape)
library(dendextend)


####PART C####

#MA------------------------- Fetching and Trimming Data Set ------------------------------------------

#MA: Date and time for when data was acquired (not provided)
Daphnia <- readDNAStringSet("SixteenS.fetch.fasta")

#MA: Converting data set into a data frame with columns containing the data label and gene sequences
dfDaphnia <- data.frame(Daphnia.Title = names(Daphnia), Daphnia.Sequence = paste(Daphnia))

#MA: Re-formatting data frame to include proper labels including species name, accession numbers(identifier), gene names and sequences
dfDaphnia$Species.Name <- word(dfDaphnia$Daphnia.Title, 2L, 3L)#MA: extracting words #2 and #3 from original label for species name
dfDaphnia$Unique.Identifier <- word(dfDaphnia$Daphnia.Title, 1L)#MA: extracting word #1 for unique identifier
dfDaphnia$Gene.Name <- word(dfDaphnia$Daphnia.Title, 6L)#MA: extracting word #6 for gene names
#dfDaphnia$gene <- str_extract(dfDaphnia$Daphnia.Title, "16S.*") MA: not necessary to include
dfDaphnia <- dfDaphnia[, c("Unique.Identifier", "Species.Name", "gene", "Gene.Name", "Daphnia.Sequence", "Daphnia.Title")]

#MA: Filtering only records that includes 16S DNA
dfDaphnia.Filter <- dfDaphnia %>%
  filter(dfDaphnia$Gene.Name == "16S")

#MA:Randomly selecting for one sample per species to reduce data size
DaphniaSub <- dfDaphnia.Filter %>%
  group_by(Species.Name) %>%
  sample_n(1)

#MA------------------------------------MSA and Phylogenetic Reconstruction----------------------------

#MA: Converting sequence data into DNAstringset for downstream analysis using function muscle from the muscle package
DaphniaSub$Daphnia.Sequence <- DNAStringSet(DaphniaSub$Daphnia.Sequence)

#MA: Multiple sequence alignment (MSA); no explanations was provided on why these parameters were used
Daphnialign <- DNAStringSet(muscle::muscle(DaphniaSub$Daphnia.Sequence, log = "log.tx", verbose = T), use.names = TRUE)

#MA: Conversion of MSA data to DNAbin to be used for calculating distance matrix
DaphniaDNAbin <- as.DNAbin(Daphnialign)
class(DaphniaDNAbin)

#MA: Calculating distance matrix using the function dist.dna from the ape package
Daphniadist <- dist.dna(DaphniaDNAbin, model = "TN93", as.matrix = TRUE, pairwise.deletion = TRUE)

#MA: Creating a distance based cluster for MSA data with the complete-linkage method; no explanation provided why parameters are chosen
Daphniacluster <- hclust(dist(Daphniadist), method = "complete")
plot(Daphniacluster)


#MA: ------------ Maverick's additions ---------------

#Creating a distance based cluster for MSA data with the single-linkage method
Daphniacluster_single <- hclust(dist(Daphniadist), method = "single")
plot(Daphniacluster)

#Creating more visually appealing phylogenetic trees

#Creating a dendrogram using 16S MSA derived distance matrix (complete-linkage)
D_complete_16S <- Daphniacluster_comp %>%
  as.dendrogram %>% #converting cluster to a dendrogram
  set("branches_k_color", k=11)  #setting color based on clusters for branches
plot(D_complete_16S)

#Creating a dendrogram using 16S MSA derived distance matrix (single-linkage)
D_single_16S <- Daphniacluster_single %>%
  as.dendrogram %>% #converting cluster to a dendrogram
  set("branches_k_color", k=11)  #setting color based on clusters for branches
plot(D_single_16S)


#Only one gene data set was used. Instead of comparing phylogenetic hypothesis of two genes, the single data set can simply be used to compare the phylogenetic hypothesis of two clustering algorithms.
#Compares single and complete linkage derived dendrograms using tanglegram from then dendextend package
tanglegram(D_complete_16S,D_single_16S,main ="Comparing Single and Complete Linkage Trees", main_left="Single-linkage", main_right = "Complete-linkage", common_subtrees_color_branches = TRUE, columns_width = c(5, 2, 5),margin_inner = 2,cex_main = .7, cex_main_left = 1, cex_main_right = 1, sort = FALSE, rank_branches = TRUE)
