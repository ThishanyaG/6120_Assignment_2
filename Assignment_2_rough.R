## 6210 Assignment 2

# The Differentiation of CytB sequences between Extinct and Extant Elephantidae species

# Load libraries 
library(rentrez)
library(BiocGenerics)
library(Biostrings)
library(tidyverse)
library(muscle)


################################################################################
# Data acquisition: get gene data using entrez search

# First, look at the databases available. A search in NCBI online shows that there is a sizable amount of information on the chosen family in the nucleotide database, so we then look at a summary of the db and the searchable parameters.
entrez_dbs()
entrez_db_summary(db = 'nucleotide')
entrez_db_searchable(db = "nucleotide")
entrez_db_links("nucleotide")

# 
elephantidae <- entrez_search(db="nucleotide", term = "Elephantidae[ORGN] AND cytochrome b")
elephantidae <- entrez_search(db="nucleotide", term = "Elephantidae[ORGN]  AND cytochrome b", retmax = elephantidae$count, use_history = T)

summary_elephant <- entrez_summary(db = "nucleotide", id = elephantidae$ids[1:300])
View(extract_from_esummary(summary_elephant, "title"))
summary_elephant1 <- entrez_summary(db = "nucleotide", id = elephantidae$ids[300:600])
View(extract_from_esummary(summary_elephant1, "title"))

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

#############################################################################################
elephant_search <- entrez_search(db = "nucleotide", term = "Elephantidae[ORGN] AND 16S AND 10:400[SLEN] NOT Bacteria[ORGN]")
elephantidae_16S<- entrez_search(db = "nucleotide", term = "Elephantidae[ORGN] AND 16S AND 10:400[SLEN] NOT Bacteria[ORGN]", retmax = elephant_search$count, use_history = T)

elephant_search1 <- entrez_search(db = "nucleotide", term = "Elephantidae[ORGN] AND cytochrome b")
elephantidae_cytb <- entrez_search(db = "nucleotide", term = "Elephantidae[ORGN] AND cytochrome b AND 500:800[SLEN]", retmax = elephant_search1$count, use_history = T)

# Remove the initial searches
rm(elephant_search, elephant_search1)


post_cytb <- entrez_post(db = "nucleotide", web_history = elephantidae_cytb$web_history)

sum_cytb1 <- entrez_summary(db = "nucleotide", elephantidae_cytb$ids[1:330])
sum_cytb2 <- entrez_summary(db = "nucleotide", elephantidae_cytb$ids[330:660])
sum_cytb3 <- entrez_summary(db = "nucleotide", elephantidae_cytb$ids[660:elephantidae_cytb$count])

View(extract_from_esummary(sum_cytb1, "title"))
View(extract_from_esummary(sum_cytb2, "title"))
View(extract_from_esummary(sum_cytb3, "title"))


set_list_16S <- c()
set_list_cytb <- c()

for (start_rec in seq(0, elephantidae_cytb$retmax, 200)) {
  fname <- paste("Elephant_cytb", start_rec, ".fasta", sep = "")
  recs <- entrez_fetch(db = "nucleotide", web_history = elephantidae_cytb$web_history, rettype = "fasta", retstart = start_rec, retmax = 200)
  write(recs, fname)
  print(paste("Wrote records to ", fname, sep = ""))
  set_list_cytb <- append(set_list_cytb, fname)
}

for (start_rec in seq(0, elephantidae_16S$retmax, 200)) {
  fname <- paste("Elephant_16S", start_rec, ".fasta", sep = "")
  recs <- entrez_fetch(db = "nucleotide", web_history = elephantidae_16S$web_history, rettype = "fasta", retstart = start_rec, retmax = 200)
  write(recs, fname)
  print(paste("Wrote records to ", fname, sep = ""))
  set_list_16S <- append(set_list_16S, fname)
}

cytb_stringSet_0 <- readDNAStringSet("Elephant_cytb0.fasta")
cytb_stringSet_200 <- readDNAStringSet("Elephant_cytb200.fasta")
cytb_stringSet_400 <- readDNAStringSet("Elephant_cytb400.fasta")
cytb_stringSet_600 <- readDNAStringSet("Elephant_cytb600.fasta")

stringSet_16S_0 <- readDNAStringSet("Elephant_16S0.fasta")
stringSet_16S_200 <- readDNAStringSet("Elephant_16S200.fasta")
stringSet_16S_400 <- readDNAStringSet("Elephant_16S400.fasta")

cytb1 <- data.frame(CytB_Title = names(cytb_stringSet_0), CytB_Sequence = paste(cytb_stringSet_0))
cytb2 <- data.frame(CytB_Title = names(cytb_stringSet_200), CytB_Sequence = paste(cytb_stringSet_200))
cytb3 <- data.frame(CytB_Title = names(cytb_stringSet_400), CytB_Sequence = paste(cytb_stringSet_400))
cytb4 <- data.frame(CytB_Title = names(cytb_stringSet_600), CytB_Sequence = paste(cytb_stringSet_600))
rm(cytb_stringSet_0, cytb_stringSet_200, cytb_stringSet_400, cytb_stringSet_600)

df_cytb1 <- merge(cytb1, cytb2, all = T)
df_cytb2 <- merge(cytb3, cytb4, all = T)
df_cytb <- merge(df_cytb1, df_cytb2, all = T)
rm(df_cytb1, df_cytb2, cytb1, cytb2, cytb3, cytb4)

e_16S1 <- data.frame(Title_16S = names(stringSet_16S_0), Sequence_16S = paste(stringSet_16S_0))
e_16S2 <- data.frame(Title_16S = names(stringSet_16S_200), Sequence_16S = paste(stringSet_16S_200))
e_16S3 <- data.frame(Title_16S = names(stringSet_16S_400), Sequence_16S = paste(stringSet_16S_400))
rm(stringSet_16S_0, stringSet_16S_200, stringSet_16S_400)

df_16S1 <- merge(e_16S1, e_16S2, all = T)
df_16S <- merge(df_16S1, e_16S3, all = T)
rm(df_16S1, e_16S1, e_16S2, e_16S3)

#############################################################################################
# Data summerization: summerise the data that has been collected




cytb_summary <- entrez_summary(db = "nucleotide", web_history = elephantidae_cytb$web_history)

?entrez_link()

ele_sum1 <- entrez_summary()
View(extract_from_esummary(ele_sum1, "title"))
ele_sum <- entrez_summary(db = "nucleotide", id = elephant_search1$ids)
View(extract_from_esummary(ele_sum, "title"))




#############################################################################################

elephas_search <-  entrez_search(db = "nucleotide", term = "Elephas[ORGN] AND cytochrome b", retmax = 100)
elephas_search
elephas_sum <- entrez_summary(db = "nucleotide", id = elephas_search$ids)
View(extract_from_esummary(elephas_sum, "title"))

loxodonta_search <-  entrez_search(db = "nucleotide", term = "Loxodonta[ORGN] AND cytochrome b", retmax = 500, use_history = T)
loxodonta_search
loxo_sum <- entrez_summary(db = "nucleotide", id = loxodonta_search$ids)
View(extract_from_esummary(loxo_sum, "title"))

mammuthus_search <-  entrez_search(db = "nucleotide", term = "Mammuthus[ORGN] AND cytochrome b", retmax = 100)
mammuthus_search
mamm_sum <- entrez_summary(db = "nucleotide", id = mammuthus_search$ids)
View(extract_from_esummary(mamm_sum, "title"))

#############################################################################################

elephant <-  entrez_search(db = "nucleotide", term = "Elephantidae[ORGN] AND 500:1000[SLEN]", retmax = 200, use_history = T)

summary <- entrez_summary(db = "nucleotide", id = elephant$ids)
el_sum <- extract_from_esummary(summary, "title")
View(el_sum)


elephas_search_16S <-  entrez_search(db = "nucleotide", term = "Elephas[ORGN] AND 16S", retmax = 100)
elephas_search_16S
elephas_sum <- entrez_summary(db = "nucleotide", id = elephas_search$ids)
View(extract_from_esummary(elephas_sum, "title"))

loxodonta_search_16S <-  entrez_search(db = "nucleotide", term = "Loxodonta[ORGN] AND 16S", retmax = 100)
loxodonta_search_16S
loxo_sum <- entrez_summary(db = "nucleotide", id = loxodonta_search$ids)
View(extract_from_esummary(loxo_sum, "title"))

mammuthus_search_16S <-  entrez_search(db = "nucleotide", term = "Mammuthus[ORGN] AND 16S", retmax = 100)
mammuthus_search_16S
mamm_sum <- entrez_summary(db = "nucleotide", id = mammuthus_search$ids)
View(extract_from_esummary(mamm_sum, "title"))


#############################################################################################

mammoth_Search <- entrez_search(db = "nucleotide", term = "Mammuthus primigenius[ORGN] AND cytb[GENE]")
mammoth_Search

mammoth_sum <- entrez_summary(db = "nucleotide", id = mammoth_Search$ids)
View(mammoth_sum)
extract_from_esummary(mammoth_sum, "title")

Asian_Search <- entrez_search(db = "nucleotide", term = "Elephas maximus[ORGN] AND cytb[GENE]")
Asian_Search
Asian_sum <- entrez_summary(db = "nucleotide", id = Asian_Search$ids)
#View(Asian_sum$)
extract_from_esummary(Asian_sum, "title")

African_Search <- entrez_search(db = "nucleotide", term = "Loxodonta Africana[ORGN] AND cytb[GENE]")
African_Search
African_sum <- entrez_summary(db = "nucleotide", id = African_Search$ids)
#View(African_sum$)
extract_from_esummary(African_sum, "title")


#############################################################################################
d_loop_search <- entrez_search(db = "nuccore", term = "Elephas maximus[ORGN] AND mitochondrial D-loop", retmax = d_loop_search$count, use_history = T)

post <- entrez_post(db = "nuccore", web_history = d_loop_search$web_history)
entrez_db_links("nuccore")

sum <- entrez_summary(db = "nucleotide", web_history = post)
View(extract_from_esummary(sum, "title"))

#d_loop_search
#dloop_fetch <- entrez_fetch(db = "nucleotide", id = d_loop_search$ids, rettype = "fasta")
#write(dloop_fetch, "dloop_fetch.fasta", sep = "\n")