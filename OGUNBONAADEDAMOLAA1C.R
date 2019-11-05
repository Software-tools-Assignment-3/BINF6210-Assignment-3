####QUESTION####
#How comparable is species richness in 2 distinct climate regions? We would be comparing 2 countries on opposite sides of the world with lot of data acquired from both locations, we would be using the vegan package as a supplement for analyzing concluding data.

####LOAD IMPORTANT LIBRARIES AND BOLD FILES####
library(tidyverse)
library(BiocManager)
library(Biostrings)
library(vegan)

CephalopodaCanada <- read_tsv("http://www.boldsystems.org/index.php/API_Public/combined?taxon=Cephalopoda&geo=Canada&format=tsv")

CephalopodaAustralia <- read_tsv("http://www.boldsystems.org/index.php/API_Public/combined?taxon=Cephalopoda&geo=Australia&format=tsv")

####ANALYZE CANADA DATA----
#We need to get an idea of the sort of data we are analyzing
nrow(CephalopodaCanada)
summary(CephalopodaCanada)
#Explore the latitude column and count the empty cell values
sum(is.na(CephalopodaCanada$lat))
#Create an histogram that shows distribution based on location
hist(CephalopodaCanada$lat)
min(CephalopodaCanada$lat, na.rm = TRUE)
#Find unique Markercode in the distribution and the most frequent marker which is COI-5P
unique(CephalopodaCanada$markercode)
which.max(table(CephalopodaCanada$markercode))
#Create subset from marker COI-5P containing nucleotide sequences
CephalopodaCanadaCOI <- CephalopodaCanada %>%
       filter(markercode == "COI-5P") %>%
       filter(str_detect(nucleotides, "[ACGT]"))
#Count how many unique BINs, marker codes and species are in the distribution
length(unique(CephalopodaCanadaCOI$bin_uri))
length(unique(CephalopodaCanadaCOI$species_name))
table(CephalopodaCanadaCOI$country)
unique(CephalopodaCanadaCOI$markercode)
#Convert the nucleotide sequences from characters to DNA sequences
class(CephalopodaCanadaCOI$nucleotides)
CephalopodaCanadaCOI$nucleotides <- DNAStringSet(CephalopodaCanadaCOI$nucleotides)
class(CephalopodaCanadaCOI$nucleotides)
#Count the occurence of each nucleotide and the proportion of the A and T nucleotides and add as a new column
CephCanNucFreq <- as.data.frame(letterFrequency(CephalopodaCanadaCOI$nucleotides, letters = c("A","C","G","T")))
CephCanNucFreq <- CephCanNucFreq %>%
  mutate(ATproportion = ((A + T) / (A + T + G + C)))
hist(CephCanNucFreq$ATproportion)
mean(CephCanNucFreq$ATproportion)
#Create a community by BIN
CephalopodaCanadaCOIcomm <- CephalopodaCanadaCOI %>%
  group_by(bin_uri) %>%
  count(bin_uri)
Canada.spread <- spread(CephalopodaCanadaCOIcomm, bin_uri, n)
#Create an accumulation curve
Can.curve <- rarecurve(Canada.spread)
#Repeat this steps on the Australian data

####ANALYZE AUSTRALIA DATA----
nrow(CephalopodaAustralia)
summary(CephalopodaAustralia)
sum(is.na(CephalopodaAustralia))
hist(CephalopodaAustralia$lat)
min(CephalopodaAustralia$lat, na.rm = TRUE)
unique(CephalopodaAustralia$markercode)
which.max(table(CephalopodaAustralia$markercode))
CephalopodaAustraliaCOI <- CephalopodaAustralia %>%
      filter(markercode == "COI-5P") %>%
      filter(str_detect(nucleotides, "[ACGT]"))
length(unique(CephalopodaAustraliaCOI$bin_uri))
length(unique(CephalopodaAustraliaCOI$species_name))
unique(CephalopodaAustraliaCOI$species_name)
unique(CephalopodaAustraliaCOI$markercode)
CephalopodaAustraliaCOI$nucleotides <- DNAStringSet(CephalopodaAustraliaCOI$nucleotides)
class(CephalopodaAustraliaCOI$nucleotides)
CephAusNucFreq <- as.data.frame(letterFrequency(CephalopodaAustraliaCOI$nucleotides, letters = c("A","C","G","T")))
CephAusNucFreq <- CephAusNucFreq %>%
  mutate(ATproportion = ((A + T) / (A + T + G + C)))
hist(CephAusNucFreq$ATproportion)
mean(CephAusNucFreq$ATproportion)
CephalopodaAustraliaCOIcomm <- CephalopodaAustraliaCOI %>%
  group_by(bin_uri) %>%
  count(bin_uri)
Australia.spread <- spread(CephalopodaAustraliaCOIcomm, bin_uri, n)
Aus.Curve <- rarecurve(Australia.spread)


####FIND COMMON SPECIES IN BOTH REGIONS----
#We use the intersect function to find a common elements in both countries' species column. As we discover, no common species native to both countries
intersect(CephalopodaCanada$species_name,CephalopodaAustralia$species_name)

#For dissimilarities between the 2 databases, we use vegdist
#First we create a new dataframe called df containing taxonomy ID column of both countries
df <- merge.data.frame(CephalopodaAustralia$order_taxID, CephalopodaCanada$order_taxID)

vegdist(df, na.rm = TRUE) #Seems too large to allocate, froze my laptop for a while

####CONCLUSION####
#There's no common species in the dataframe.
#The accumulation curve is close 