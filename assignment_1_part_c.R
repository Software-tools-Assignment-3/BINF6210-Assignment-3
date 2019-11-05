#Part C from assignment 1
#Importing data
chytridiomycota <- read_tsv('Chytridiomycota.txt')
#Lets first explor the data by looking at the columns
names(chytridiomycota)
#Now lets see how many countries were specimens collected from
length(unique(chytridiomycota$country))
#Which country do you think collected the most specimens
names(which.max(table(chytridiomycota$country)))
#To begin with lets look at how many unique order, family, genus and species have been collected
length(unique(chytridiomycota$order_name))
length(unique(chytridiomycota$family_name))
length(unique(chytridiomycota$genus_name))
length(unique(chytridiomycota$species_name))
#What is the least frequent species
names(which.min(table(chytridiomycota$species_name)))
#Now we are going to determine the habitat of the most frequent species
names(which.max(table(chytridiomycota$species_name)))
subset(chytridiomycota, chytridiomycota$species_name == 'Batrachochytrium dendrobatidis', select = chytridiomycota$habitat)
#Now lets make sure each barcode sequence has a unique sequenceID
length(chytridiomycota$nucleotides) == length(chytridiomycota$sequenceID)
#Now I am going further and draw rarefaction curves
#I am going to group the data by the number of records in each country
chytridiomycota_count_country <- chytridiomycota %>%
  group_by(country) %>%
  count(country)
#Now I will use spread() function to create a community object
chytridiomycota_spread <- spread(chytridiomycota_count_country, country, n)
#Now lets plot it
chytridiomycota_accum <- rarecurve(chytridiomycota_spread)
#Another exploration would be to answer: How many unique species are there in the most frequenct genus
#In order to answer this question I am going to use a function called filter()
most_fre_genus <- names(which.max(table(chytridiomycota$genus_name)))
filt_chytridio <- filter(chytridiomycota, genus_name == most_fre_genus)
filt_chytridio$species_name
#Now lets calculate the mean length of each barcode sequence belonging to the most and the least frequent species
nchar(na.omit(chytridiomycota$nucleotides[chytridiomycota$species_name == 'Alogomyces tanneri']))
nchar(na.omit(chytridiomycota$nucleotides[chytridiomycota$species_name == 'Batrachochytrium dendrobatidis']))
#Now lets calculate the count of substring GC of the most frequent species
chytridiomycota$nucleotides <- DNAStringSet(chytridiomycota$nucleotides)
letterFrequency(chytridiomycota$nucleotides, letters = 'GC')
