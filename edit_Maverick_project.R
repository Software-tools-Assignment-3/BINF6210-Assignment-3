#This script has been edited by Pasha Talebi Charmchi

#PART C - Adapting Concepts and Skills to Create New Research 

#Project topic and objectives:

#In this mini-project, I wanted to determine how much of impact different countries has contributed to the efforts of classifying all the members of the Sipuncula taxa. For this, I need to know which countries have submitted barcodes and how many have they submitted at the current time. To get the bigger picture, I hope to determine the current status of the classification globally and then follow up with how each country have contributed to the effort. In this respect, I want to determine the genetic diversity of species found in different geographic locations. Ultimately, I hope to determine which country had the highest impact based on the genetic diversity of the barcodes recorded and submitted from said country.


################################################################################################################

#Package initiations
library(tidyverse)
library(vegan)
library(Biostrings)

################################################################################################################

#Taxon of Sipuncula, files acquired Oct 2, 2019 at 11:30 AM.
Sipuncula <- read_tsv("http://www.boldsystems.org/index.php/API_Public/combined?taxon=Sipuncula&format=tsv")

write_tsv(Sipuncula, "Sipuncula.tsv")

################################################################################################################

#Checking basic attributes of Sipincula data set
class(Sipuncula)
summary(Sipuncula)
names(Sipuncula)

################################################################################################################

#Checking for geographic distributions

#Checking missing latitudal records
sum(is.na(Sipuncula$lat))
hist(Sipuncula$lat, xlab = 'Latitude', ylab = 'Records', main = 'Fig 1: Distribution of Records by Latitude')
#My peer used above function to check geographic distributions which is actually a built-in function for data visualization
#My alternative way to plot a histogram would be to use ggplot as the output would be more beautiful than using his() function. Despite its complexity we can easily add more features to our plot by adding more layers.
ggplot(data = Sipuncula, aes(x = Sipuncula$lat)) + geom_histogram(color = 'darkgray', fill = 'white', binwidth = 60) + labs(x = 'Latitude', y = 'Records', title = 'Distribution of Records by Latitude')

################################################################################################################

#Checking basic parameters
mean(Sipuncula$lat, na.rm = TRUE)
min(Sipuncula$lat, na.rm = TRUE)
max(Sipuncula$lat, na.rm = TRUE)
#My peer used above functions to check basic attributes 
#My alternative way would be to use summary() function as it is more efficient and we can achieve the same results by just writing one line of code.
summary(Sipuncula$lat)

################################################################################################################

#Checking how many records are there for each country

sum(is.na(Sipuncula$country)) # Quick check for records with no entries for country of origin

Sipuncula.per.country.count <- Sipuncula %>%
  #Groups the records by country
  group_by(country) %>%
  #Counts the number of records within each country and returns the counts in a descending order
  summarize(no.records = length(processid)) %>%
  arrange(desc(no.records)) %>%
  print()

################################################################################################################

#Filtering data set

#Determing which genetic markers was used most frequently
Sipuncula.markercode <- Sipuncula %>%
  #All records with no nucleotide data will be filtered out.
  filter(!is.na(nucleotides)) %>%
  #Group data by their genetic markers.
  group_by(markercode) %>%
  #Summarizes data for each genetic marker and length() counts how many records are there per genetic marker.
  summarize(n = length(processid)) %>%
  #Return the counts for each genetic marker in descending order.
  arrange(desc(n)) %>%
  print()
###############################################################################################################

#It would be wise to use piping but an alternative way would be to use table() function. By using table() function we can achieve the same result with just one line of code.
table(Spiuncula$markercode)

###############################################################################################################

#Creating a data subset consisting only of records with nucleotide sequences from the genetic marker COI-5P. Only this genetic marker has relevant number of records. 
Sipuncula.COI_5P <- Sipuncula %>%
  #Removes all records without nucleotide data
  filter(!is.na(nucleotides)) %>%
  #Takes all records that corresponds to COI-5P genetic marker
  filter(markercode == "COI-5P") %>%
  #Takes all records that contains the four DNA bases
  filter(str_detect(nucleotides, "[ACGT]"))

################################################################################################################

#Creating a species accumulation curves for Sipuncula records with nucleotide sequences from the genetic marker COI-5P. 

#Global accumulation curve: this curve will demonstrate whether classification of all species of the Taxon Sipuncula is near completion.

#Groups all records by BINs and put them all in one global site and counts how many records exist for each BIN.
Sipuncula.COI_5P.global <- Sipuncula.COI_5P %>%
  group_by(bin_uri) %>%
  count(bin_uri)

#Converts commAssign1 into a community object. The BINs will become column headers and the record counts for each bin will be the row values (the row just being the one global site).
Sipuncula.COI_5P.commglobal <- spread(Sipuncula.COI_5P.global, bin_uri, n)

#Will create an accumulation curve (using function rarecurve()) that treats all BINs as a function of one site, global. 
Sipuncula.globalmodel <- rarecurve(Sipuncula.COI_5P.commglobal, xlab = 'No. of Records', ylab = 'No. of Species (BINs)', main ='Fig 2: Sipuncula Species Accumulation Curve in a Single Global Site')

#Species accumulation curve using BINs as functions of their countries of discovery. This curve may demonstrate patterns in rate of discovery as more countries (or sites) get added in. Can answer whether each country adds alot of unique species or there is a tendency of different countries having the same species. 

#Formatting data to create a data frame that returns the number of samples per BIN in each country.

#Will group data by BINs by country to allow the count for the number of specimens that have each BINs per country
Sipuncula.COI_5P.country <- Sipuncula.COI_5P %>%
  group_by(country, bin_uri) %>% 
  count(bin_uri)

#Remove all records that doesn't have entries for their country of origin and BINs. The goal is to see the effects of addition of countries to the number of unique BINs. A missing record of either makes the data irrelevant.
Sipuncula.COI_5P.country.rm <- Sipuncula.COI_5P.country %>%
  filter(!is.na(country)) %>%
  filter(!is.na(bin_uri))

#Converting data frame containing the sample counts per BINs per country into a community object using spread()
Sipuncula.COI_5P.comm.country <- spread(Sipuncula.COI_5P.country.rm, bin_uri, n)

#Converts NAs to zeroes, to tell R that barcoding a particular BIN have yet to be done in that country.Required for defining the ylim of the plot.
Sipuncula.COI_5P.comm.country[is.na(Sipuncula.COI_5P.comm.country)] <- 0

#Sets country as a row rather than a data column (to use as the X-axis in a species accumulation curve).
Sipuncula.COI_5P.comm.country <- Sipuncula.COI_5P.comm.country %>%
  remove_rownames %>%
  column_to_rownames(var="country")

#Creates a species accumulation curve with numbers of species as the function of the number of countries.
Sipuncula.COI_5P.countrymodel <- specaccum(Sipuncula.COI_5P.comm.country)

#Plots the model curve
plot(Sipuncula.COI_5P.countrymodel, xlab = 'No. of Sites (Countries)', ylab = 'No. of Species (BINs)', main ='Fig 3: Sipuncula Species Accumulation Curve from Multiple Sites ')

#Dataframe attributes check
summary(Sipuncula.COI_5P.comm.country)
names(Sipuncula.COI_5P.comm.country)
str(Sipuncula.COI_5P.comm.country)

###################################################################################################################################

#Creating cluster to determine dissimilarity between the species from each country.

#Calculates relative distances among samples using Bray-Curtis method
Sipuncula.bray <- vegdist(Sipuncula.COI_5P.comm.country, method = "bray", rm = TRUE)

#Uses average-linkage algorithm to cluster communities
Sipuncula.brayclust <- hclust(Sipuncula.bray, method = "average")

#returns a community cluster diagram
plot(Sipuncula.brayclust, ylab = "Bray-Curtis Dissimilarity", main = '', xlab = "Fig 4: Genomic dissimilarities between species from different countries")


##############################################################################################################

#Determining diversity of species within countries with the most number of records.

#Creates a subset that contains counts of specimens per BINs per country
Sipuncula.COI_5P.Top.5 <- Sipuncula %>%
  filter(country == "Russia" | country == "China"| country == "United States"| country == "Canada"| country == "Thailand") %>% #extracts records for top 5 countries with highest number of records
  group_by(country, bin_uri) %>% #groups by BINs per country
  count(bin_uri)

#Removes all missing values in variables country and bin_uri
Sipuncula.COI_5P.Top.5.rm <- Sipuncula.COI_5P.Top.5 %>%
  filter(!is.na(country)) %>%
  filter(!is.na(bin_uri))

#Transforms data frame into a community object
Sipuncula.COI_5P.comm.Top.5 <- spread(Sipuncula.COI_5P.Top.5.rm, bin_uri, n)

#Converts NAs to zeroes for downstream processing
Sipuncula.COI_5P.comm.Top.5[is.na(Sipuncula.COI_5P.comm.Top.5)] <- 0

#Sets country as a row rather than a data column 
Sipuncula.COI_5P.comm.Top.5 <- Sipuncula.COI_5P.comm.Top.5 %>%
  remove_rownames %>%
  column_to_rownames(var="country")

#returns the Simpson's Reciprocal Index for each country. Will show which country has the highest diversity
diversity(Sipuncula.COI_5P.comm.Top.5, index = "invsimpson")


#returns the number of species 
specnumber(Sipuncula.COI_5P.comm.Top.5)

##########################################################################################################

#There were 21 different countries that have submitted records, however most of the BINs did not have entries for their countries of origins. The sampling of Sipuncula members globally is largely incomplete (Fig 2.) and it is apparent that there are still several unknown species within each country listed in this study and those which are not (Fig 3.). The species found in North America, Japan and French Polynesia have been observed to be more genetically different relative to the rest of the species found in other locations (Fig 4.). Most of the records come from Russia, however these records belonged only to 3 unique species, suggesting a low genetic diversity. The country that have the highest genetic diversity is the United States (17 unique species), which agrees with the finding of species from North America being more genetically different. In conclusion, sampling of Sipuncula is largely incomplete, and many species are yet to be sampled within most countries. 

##########################################################################################################

#As a whole, the script is well commented and organized. However, providing some information about the reasons for using the function used in this project would improve the script and it would be better if the comments were put in front of each code rather than in a new line. 
