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

####PART C####
Daphnia <- readDNAStringSet("SixteenS.fetch.fasta")

dfDaphnia <- data.frame(Daphnia.Title = names(Daphnia), Daphnia.Sequence = paste(Daphnia))

dfDaphnia$Species.Name <- word(dfDaphnia$Daphnia.Title, 2L, 3L)
dfDaphnia$Unique.Identifier <- word(dfDaphnia$Daphnia.Title, 1L)
dfDaphnia$Gene.Name <- word(dfDaphnia$Daphnia.Title, 6L)
dfDaphnia$gene <- str_extract(dfDaphnia$Daphnia.Title, "16S.*")
dfDaphnia <- dfDaphnia[, c("Unique.Identifier", "Species.Name", "gene", "Gene.Name", "Daphnia.Sequence", "Daphnia.Title")]
dfDaphnia.Filter <- dfDaphnia %>%
  filter(dfDaphnia$Gene.Name == "16S")

DaphniaSub <- dfDaphnia.Filter %>%
  group_by(Species.Name) %>%
  sample_n(1)
DaphniaSub$Daphnia.Sequence <- DNAStringSet(DaphniaSub$Daphnia.Sequence)
Daphnialign <- DNAStringSet(muscle::muscle(DaphniaSub$Daphnia.Sequence, log = "log.tx", verbose = T), use.names = TRUE)
DaphniaDNAbin <- as.DNAbin(Daphnialign)
class(DaphniaDNAbin)
Daphniadist <- dist.dna(DaphniaDNAbin, model = "TN93", as.matrix = TRUE, pairwise.deletion = TRUE)
Daphniacluster <- hclust(dist(Daphniadist), method = "complete")
plot(Daphniacluster)