#This R script has been edited by Pasha Talebi Charmchi

# Load Packages
library(tidyverse)
library(muscle)
library(Biostrings)
library(DECIPHER)
library(phytools)
library(ape)
library(phangorn)
library(dendextend)

#Order of Equisetales , files acquired Oct 30, 2019 at 3:00 PM.
Equisetales <- read_tsv("http://www.boldsystems.org/index.php/API_Public/combined?taxon=Equisetales&format=tsv")

#Writing data set into a file
write_tsv(Equisetales, "Equisetales.tsv")

Equisetales <- read_tsv("Equisetales.tsv")

################################################################################################################

##Data attributes and cleaning data

#Checking basic attributes of Sipincula data set
class(Equisetales)
summary(Equisetales)
names(Equisetales)

#extracting specific columns useful for phylogenetic studies
Equisetales_filtered <- Equisetales %>%
  select(processid, bin_uri, genus_name, species_name, country, lat, lon, markercode, nucleotides)

#Checking how many unique species in data set
length(unique(Equisetales_filtered$species_name))

#Checking how many unique genus in data set
length(unique(Equisetales_filtered$genus_name))

#Quick check for records with no entries for country of origin
sum(is.na(Equisetales_filtered$country)) 
#Quick check for records with no entries for latitude
sum(is.na(Equisetales_filtered$lat))
#Quick check for records with no entries for longitude
sum(is.na(Equisetales_filtered$lon))

#Checking for how many countries have records of species
length(unique(Equisetales_filtered$country))

#Counts the number of records within each country and returns the counts in a descending order
Equ_per_country <- Equisetales_filtered %>%
  group_by(country) %>%
  summarize(no.records = length(processid)) %>%
  arrange(desc(no.records)) %>%
  print()

#Removing records with missing data for species name, genus name, markercode, country, latitude, longitude and missing nucleotides
Clean_Equ <- Equisetales_filtered %>%
  filter(!is.na(species_name)) %>%
  filter(!is.na(genus_name)) %>%
  filter(!is.na(country)) %>%
  filter(!is.na(lat)) %>%
  filter(!is.na(lon)) %>%
  filter(!is.na(markercode)) %>%
  filter(!is.na(nucleotides)) 

#Attribute check for cleaned up dataframe
dim(Clean_Equ)
names(Clean_Equ)
unique(Clean_Equ$species_name)
sum(is.na(Clean_Equ$nucleotides))


################################################################################################################

#Checking and filtering nucleotide data

#Checking basic attributes
min(str_length(Clean_Equ$nucleotides))
mean(str_length(Clean_Equ$nucleotides))
median(str_length(Clean_Equ$nucleotides))
max(str_length(Clean_Equ$nucleotides))
#My peer used above functions to check basic attributes for nucleotide sequences he retrieved from data base
#My alternative way to check basic attributes for nucleotides data would be to use both nchar() and summary() function as it is more efficient and we can achieve the same results by just writing one line of code.
summary(nchar(Clean_Equ$nucleotides))


#Checking distribution of sequence lengths
hist(str_length(Clean_Equ$nucleotides))
#My peer used above function to check the distribution of sequence length which is actually a basic function for data visualization
#My alternative way to plot a histogram would be to use ggplot2 as the output would be more beautiful than a basic built-in function. Despite its complexity we can easily add features to our plot by adding more layers.
ggplot(data = Clean_Equ, aes(x = nchar(Clean_Equ$nucleotides))) + geom_histogram(binwidth = 60) + labs(x = 'Sequences', y = 'Number Of Nucleotides', title = 'Distribution Of Sequence Lengths')

#Histogram displays stragglers in terms of sequence length; removing sequences < 550 bp and sequences > 600
#Removing sequences of these lengths will help in removing outliers

#Removing sequences greater than 600 bp
Equ_trimmed <- Clean_Equ %>%
  filter(!(str_length(Clean_Equ$nucleotides) > 600)) 

#Removing sequences less than 550 bp 
Equ_trimmed <- Equ_trimmed %>%
  filter(!(str_length(Equ_trimmed$nucleotides) < 550)) 

#Removing an outlier found from clustering, and MSA was re-ran
Equ_trimmed <- Equ_trimmed[-12,]

#Checking number of unique species
length(unique(Equ_trimmed$species_name))

#Checking basic attributes
min(str_length(Equ_trimmed$nucleotides))
mean(str_length(Equ_trimmed$nucleotides))
median(str_length(Equ_trimmed$nucleotides))
max(str_length(Equ_trimmed$nucleotides))
hist(str_length(Equ_trimmed$nucleotides))

################################################################################################################

#Reconstructing Phylogeny

#Reformatting sequences for downstream analysis with muscle package
Equ_trimmed$nucleotides <- DNAStringSet(Equ_trimmed$nucleotides)

#Setting processid as identifiers for each sequence
names(Equ_trimmed$nucleotides) <- Equ_trimmed$species_name

#Check if reformatting was done correctly
class(Equ_trimmed$nucleotides)

#------------------------Multiple Sequence Alignment (MSA)---------------------------------

#Using muscle for MSA with default parameters
Equ.alignment <- DNAStringSet(muscle::muscle(Equ_trimmed$nucleotides), use.names = TRUE)

#Checking Alignment
BrowseSeqs(Equ.alignment)

#calculating mean, max, min and observing distribution of gaps through a histogram
min(unlist(lapply(Equ.alignment, str_count, "-")))
max(unlist(lapply(Equ.alignment, str_count, "-")))
mean(unlist(lapply(Equ.alignment, str_count, "-")))
hist(unlist(lapply(Equ.alignment, str_count, "-")))

#Converting data type into DNAbin for more downstream analysis
dnaBin.Equ <- as.DNAbin(Equ.alignment)

#Creating a distance matrix using the "K80" model
DM.Equ <- dist.dna(dnaBin.Equ, model = "K80", as.matrix = TRUE, pairwise.deletion = TRUE)
#My peer used below functions to compute a distance matrix needed for clustering sequences
#My alternative way to create a distance matrix would be to use kdistance() function
library(kmer)
kdis <- kdistance(x = as.DNAbin(Equ.alignment), k = 3)


#Clustering using the UPGMA method with a 2% divergence threshold 
clusters.Equ <- IdClusters(DM.Equ,
                           method = "UPGMA",
                           cutoff= 0.02,      #2% sequence difference as species boundary
                           showPlot = TRUE,
                           type = "clusters", #returns a clusters output
                           verbose = TRUE)
#My alternative way to cluster our sequences would be as follow
hclust(kdis, "average")

#Reconstructing phylogenetic relationship using Neighbourhood-joining method from ape package
Equ_NJ <- nj(DM.Equ) #NJ from ape package 

#creating a phylogentic tree with Neibough Joining data
plotTree(Equ_NJ, fsize = 0.6, pts = TRUE)

##########################################################################################################

#I will be taking one sample per species to create a better (simpler) visualization of relationships within the organisms within Equisetales.


#Combining original filtered search results with their cluster affiliations
Equ_OG_clusters <- cbind(Equ_trimmed, clusters.Equ[[1]])

#Converting data type of cluster from integer into factor
Equ_OG_clusters$cluster <- as.factor(Equ_OG_clusters$cluster)

#Data check to see if conversion properly applied
class(Equ_OG_clusters$cluster)

#Noticed that two of the same species were being included when I wanted to take only one record per species (for next step).This command changes the name to match the "Genus species" format.
Equ_OG_clusters$species_name <- str_replace(Equ_OG_clusters$species_name, "subsp\\.* affine", "")


#Random selection of one representative sequence per cluster
Equ_SeqPer_Species <- Equ_OG_clusters %>%
  group_by(species_name) %>%
  sample_n(1)


#checking the number of unique species
length(unique(Equ_SeqPer_Species$species_name))

#checking the number of clusters
length(unique(Equ_SeqPer_Species$cluster))

#converting sequences into a DNAStringSet data type
Equ_SeqPer_Species$nucleotides <- DNAStringSet(Equ_SeqPer_Species$nucleotides)

#Setting species name as identifiers for each sequence
names(Equ_SeqPer_Species$nucleotides) <- Equ_SeqPer_Species$species_name


#Performing multiple sequence alignment of each representative sequence from each species
Equ.per.species.alignment <- DNAStringSet(muscle::muscle(Equ_SeqPer_Species$nucleotides), use.names = TRUE)


#converting alignment data into a DNAbin
Per.species.DNAbin <- as.DNAbin(Equ.per.species.alignment)

#Creating a distance matrix using the K80 model
DM.per.species <- dist.dna(Per.species.DNAbin, model = "K80", as.matrix = TRUE, pairwise.deletion = TRUE)


#I wanted to compare the phylogentic trees from single and and average-linkage clustering. I chose to compare these two particular clustering algorithms because average-linkage aims to strike a balance between single and complete-linkage clustering. Doing a mirror comparison will allow me to determine how different or similar they are in producing a phylogenetic tree for the same data set. 


#Creating a single-linkage clustering based tree using the upgma function from the phangorn package
SUPGMAEq <- upgma(DM.per.species, method = "single")
plotTree(SUPGMAEq)
#Labelling the tips with the clusters where each species belong to
tiplabels(Equ_SeqPer_Species$`clusters.Equ[[1]]`,adj=c(0.9,0.5))

#Creating an average-linkage based tree using the upgma function from the phangorn package
CUPGMAEq <- upgma(DM.per.species, method = "average")
plotTree(CUPGMAEq)
#Labelling the tips with the clusters where each species belong to
tiplabels(Equ_SeqPer_Species$`clusters.Equ[[1]]`,adj=c(0.9,0.5))


#Comparing the dendograms based on single-linkage and average-linkage clustering.
tanglegram(SUPGMAEq,CUPGMAEq, main_left="Single", main_right = "Average",common_subtrees_color_branches = TRUE, columns_width = c(4, 3, 4),margin_inner = 10)

#As a whole, the script is well commented. However, providing some information about the puporse of the project and the reasons for using the function used in this project would improve the script
