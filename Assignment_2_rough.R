## 6210 Assignment 2

# The Differentiation of CytB sequences between Extinct and Extant Elephantidae species

# Load libraries 
library(rentrez)
library(BiocGenerics)
library(Biostrings)
library(tidyverse)
library(muscle)
library(DECIPHER)


################################################################################
# Data acquisition: get gene data using entrez search

# First, look at the databases available. A search in NCBI online shows that there is a sizable amount of information on the chosen family in the nucleotide database, so we then look at a summary of the db and the searchable parameters.
entrez_dbs()
entrez_db_summary(db = 'nucleotide')
entrez_db_searchable(db = "nucleotide")
entrez_db_links("nucleotide")

# Do a general search and summary of elephantidae cytochrome B data available. 
elephantidae <- entrez_search(db="nucleotide", term = "Elephantidae[ORGN] AND cytochrome b")
elephantidae <- entrez_search(db="nucleotide", term = "Elephantidae[ORGN]  AND cytochrome b", retmax = elephantidae$count, use_history = T)

summary_elephant <- entrez_summary(db = "nucleotide", id = elephantidae$ids[1:300])
View(extract_from_esummary(summary_elephant, "title"))
summary_elephant1 <- entrez_summary(db = "nucleotide", id = elephantidae$ids[300:600])
View(extract_from_esummary(summary_elephant1, "title"))

# Remove objects
rm(elephantidae, summary_elephant, summary_elephant1)

# Do a search for the chosen family at a given sequence length. determine the number of data points found in the initial search and use that to do another search that includes all the sequences available.

# Extinct species:
extinct_elephants <- entrez_search(db = "nucleotide", term = "Elephas antiquus[ORGN] OR Mammut[ORGN] OR Mammuthus[ORGN] OR Elephas cypriotes[ORGN] OR Elephas maximus asurus[ORGN] AND cytochrome b AND 400:1000[SLEN]")

extinct_elephants <- entrez_search(db = "nucleotide", term = "Elephas antiquus[ORGN] OR Mammut[ORGN] OR Mammuthus[ORGN] OR Elephas cypriotes[ORGN] OR Elephas maximus asurus[ORGN] AND cytochrome b AND 400:1000[SLEN]", retmax = extinct_elephants$count, use_history = T)

extinct_summary <- entrez_summary(db = "nucleotide", id = extinct_elephants$ids[1:300])
View(extract_from_esummary(extinct_summary, "title"))

# Extant species:
extant_elephants <- entrez_search(db = "nucleotide", term = "Loxodonta africana[ORGN] OR Loxodonta cyclotis OR (Elephas maximus[ORGN] NOT Elephas maximus asurus[ORGN]) AND cytochrome b AND 500:1000[SLEN]")

extant_elephants <- entrez_search(db = "nucleotide", term = "Loxodonta africana[ORGN] OR Loxodonta cyclotis OR (Elephas maximus[ORGN] NOT Elephas maximus asurus[ORGN]) AND cytochrome b AND 500:1000[SLEN]", retmax = extant_elephants$count, use_history = T)

extant_summary <- entrez_summary(db = "nucleotide", id = extant_elephants$ids[1:300])
View(extract_from_esummary(extant_summary, "title"))

# Remove summaries
rm(extinct_summary, extant_summary)
################################################################################
#Fetch data: get the fasta files of the data and create data frames of the sequences

# Fetch data and create string sets of the fasta files 
recs_extinct <- entrez_fetch(db = "nucleotide", web_history = extinct_elephants$web_history, rettype = "fasta", retmax = extinct_elephants$count)
write(recs_extinct, "extinct_elephant_cytb.fasta")
extinct_cytb_stringSet <- readDNAStringSet("extinct_elephant_cytb.fasta")

# Make data frames of the string sets
df_extinct <- data.frame(Title = names(extinct_cytb_stringSet), Sequence = paste(extinct_cytb_stringSet))

# Add another column to the data frame that contains the species name of the data point. 
df_extinct$Species_Name <- word(df_extinct$Title, 2L, 3L)
df_extinct <- df_extinct[, c("Title", "Species_Name", "Sequence")]
View(df_extinct)


# And repeat for extant species
recs_extant <- entrez_fetch(db = "nucleotide", web_history = extant_elephants$web_history, rettype = "fasta", retmax = extant_elephants$count)
write(recs_extant, "extant_elephant_cytb.fasta")
extant_cytb_stringSet <- readDNAStringSet("extant_elephant_cytb.fasta")

df_extant <- data.frame(Title = names(extant_cytb_stringSet), Sequence = paste(extant_cytb_stringSet))

# Add another column to the data frame that contains the species name of the data point. 
df_extant$Species_Name <- word(df_extant$Title, 2L, 3L)
df_extant <- df_extant[, c("Title", "Species_Name", "Sequence")]
View(df_extant)

################################################################################
# Summarize the data
summary(df_extinct)
summary(df_extant)

# check to make sure there are no NAs in the sequence field for both data frames
df_extinct %>%
  count(is.na(Sequence))
df_extant %>%
  count(is.na(Sequence))

# count the different instances of each species in both data frames
df_extinct %>%
  group_by(Species_Name) %>%
  count(Species_Name)

df_extant %>%
  group_by(Species_Name) %>%
  count(Species_Name)

# look at the summary for the length of sequences in each data frame
summary(nchar(df_extinct$Sequence))
summary(nchar(df_extant$Sequence))

# create a histogram of the sequence lengths in each data frame
hist(nchar(df_extinct$Sequence), xlab = "Sequence Length", ylab = "
     Frequency", main = "Frequency Histogram of CytB Sequence Lengths in Extinct Species")

hist(nchar(df_extant$Sequence), xlab = "Sequence Length", ylab = "
     Frequency", main = "Frequency Histogram of CytB Sequence Lengths in Extant Species")

# Looking at the length of all the sequences in each set
df_extinct %>%
  count(nchar(Sequence))

df_extant %>%
  count(nchar(Sequence))

################################################################################
# Alignment
?muscle

extinct_alignment <- DNAStringSet(muscle::muscle(extinct_cytb_stringSet))
extant_alignment <- DNAStringSet(muscle::muscle(extant_cytb_stringSet))

BrowseSeqs(extinct_alignment)
BrowseSeqs(extant_alignment)

################################################################################



