#Part C
####QUESTION----
#I aimed to construct a pylogenetic relationship among Staphylococcus species. I first retrieved sequences from NCBI. I used 16S gene as it is the most common method used for Staphylococcus identification Then I create a datafram with column names that includes my sequences and their titles. After checking my sequences I aligned them using Muscle function so that I can compute a distance matrix based on k-mers and finally using hclust function I clustered them and drew a plot.
#####LOAD NECESSARY PACKAGES----
library(tidyverse)
library(stringi)
library(ape)
library(Biostrings)
library(muscle)
library(msa)
library(DECIPHER)
library(rentrez)
library(seqinr)
library(kmer)
library(dendextend)

#####SOLUTION PLUS EDITS----
#Search and download Staphylococcus, 16S genes from the NCBI nucleotide database
staphylococcus_search1 <- entrez_search(db = "nuccore", term = "(staphylococcus[ORGN] AND 16S[TITL]) NOT (genome[TITL])", retmax = 10)
staphylococcus_fetch <- entrez_fetch(db = "nuccore", id = staphylococcus_search1$ids, rettype = "fasta")

#Write data on disc in fasta format
write(staphylococcus_fetch, "staphylococcus_fetch.fasta", sep = "\n")

#We need to convert the data to DNA sequences 
staphylococcus_stringSet <- readDNAStringSet("staphylococcus_fetch.fasta")

#Convert the collected data into a dataframe for easy manipulation in R
dfstaphylococcus <- data.frame(staphylococcus_Title = names(staphylococcus_stringSet), staphylococcus_Sequence = paste(staphylococcus_stringSet)) 

#Some simple explorations
class(dfstaphylococcus$staphylococcus_Sequence)
is.na(dfstaphylococcus$staphylococcus_Sequence)
summary(nchar(as.character(dfstaphylococcus$staphylococcus_Sequence)))

#Create a more descriptive data set by adding more columns
dfstaphylococcus$Species_Name <- word(dfstaphylococcus$staphylococcus_Title, 2L, 3L)
dfstaphylococcus$Unique_identifier <- word(dfstaphylSococcus$staphylococcus_Title, 1L)
dfstaphylococcus$Gene_name <- word(dfstaphylococcus$staphylococcus_Title, 6L)
dfstaphylococcus$gene <- str_extract(dfstaphylococcus$staphylococcus_Title, "16S.*")
dfstaphylococcus <- dfstaphylococcus[, c("Unique_identifier", "Species_Name", "gene", "Gene_name", "staphylococcus_Sequence", "staphylococcus_Title")]

#more exploration, find the number of total sequences and unique species and create a histogram
length(dfstaphylococcus$staphylococcus_Sequence)
length(unique(dfstaphylococcus$Species_Name))
hist(str_count(dfstaphylococcus$staphylococcus_Sequence))

#Perform a multiple sequence analysis, we use kmer for determining distance
aln_staph <- muscle::muscle(DNAStringSet(dfstaphylococcus$staphylococcus_Sequence))
staphylococcus_dis <- kdistance(x = as.DNAbin(aln_staph), k = 3)

#Create distance based cluster from alignment data using single linkage
staph_single <- hclust(dist(staphylococcus_dis), method = "single")
plot(staph_single)

#Create distance based cluster from alignment data using complete linkage
staph_complete <- hclust(dist(staphylococcus_dis), method = "complete")
plot(staph_complete)

#convert the clusters to dendrogram, single linkage
staph16Ssing <- staph_single %>%
  as.dendrogram %>% 
  set("branches_k_color", k=10)
plot(staph16Ssing)

#Convert clusters to dendrogram, complete linkage
staph16Scomp <- staph_complete %>%
  as.dendrogram %>% #converting cluster to a dendrogram
  set("branches_k_color", k=10)
plot(staph16Scomp)

#Create a tanglegram from the complete and single linkage using the dendextend package
tanglegram(staph16Scomp, staph16Ssing, "Comparison of Single and complete linkage trees", main_left = "Single", main_right = "Complete", common_subtrees_color_branches = T, columns_width = c(5, 3, 5), margin_inner = 2, cex_main = .6, cex_main_left = 1, cex_main_right = 1, sort = F, rank_branches = T)

#Use browseseq to view alignment
BrowseSeqs(DNAStringSet(aln_staph))

####CONCLUSION----
#one of my sequences seemed to be distantly related to the other sequences so I decided to Blast it.It was confirmed that it belongs to Staphylococcus aureus strain FC2929 so I kept it. As a whole, my sequences were divided into two major groups and one sequences with a distant relationship. However, I would love to compare my clusters with those determined by DNA-DNA hybridization to see if they are consistent.
