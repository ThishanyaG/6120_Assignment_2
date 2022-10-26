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

elephantidae <- entrez_search(db="nucleotide", term = "Elephantidae[ORGN]")
elephantidae <- entrez_search(db="nucleotide", term = "Elephantidae[ORGN]", retmax = elephantidae$count, use_history = T)

summary_elephant <- entrez_summary(db = "nucleotide", id = elephantidae$ids[1:300])
View(extract_from_esummary(summary_elephant, "title"))

summary_elephant1 <- entrez_summary(db = "nucleotide", id = elephantidae$ids[300:600])
View(extract_from_esummary(summary_elephant1, "title"))


elephantidae <- entrez_search(db="nucleotide", term = "Elephantidae[ORGN] NOT PREDICTED[WORD]")
elephantidae <- entrez_search(db="nucleotide", term = "Elephantidae[ORGN]", retmax = elephantidae$count, use_history = T)
elephantidae$count

summary_elephant <- entrez_summary(db = "nucleotide", id = elephantidae$ids[1:300])
View(extract_from_esummary(summary_elephant, "title"))


# Do a search for the chosen family at a given sequence length. determine the number of data points found in the initial search and use that to do another search that includes all the sequences available. This was done for the genes: 16S and cytochrome B.

elephant_search <- entrez_search(db = "nucleotide", term = "Elephantidae[ORGN] AND 16S AND 10:400[SLEN] NOT Bacteria[ORGN]")
elephantidae_16S<- entrez_search(db = "nucleotide", term = "Elephantidae[ORGN] AND 16S AND 10:400[SLEN] NOT Bacteria[ORGN]", retmax = elephant_search$count, use_history = T)

elephant_search1 <- entrez_search(db = "nucleotide", term = "Elephantidae[ORGN] AND cytochrome b AND 500:800[SLEN]")
elephantidae_cytb <- entrez_search(db = "nucleotide", term = "Elephantidae[ORGN] AND cytochrome b AND 500:800[SLEN]", retmax = elephant_search1$count, use_history = T)

# Remove the initial searches
rm(elephant_search, elephant_search1)

#######################################################################################
post_cytb <- entrez_post(db = "nucleotide", web_history = elephantidae_cytb$web_history)

sum_cytb1 <- entrez_summary(db = "nucleotide", elephantidae_cytb$ids[1:330])
sum_cytb2 <- entrez_summary(db = "nucleotide", elephantidae_cytb$ids[330:660])
sum_cytb3 <- entrez_summary(db = "nucleotide", elephantidae_cytb$ids[660:elephantidae_cytb$count])

View(extract_from_esummary(sum_cytb1, "title"))
View(extract_from_esummary(sum_cytb2, "title"))
View(extract_from_esummary(sum_cytb3, "title"))

################################################################################
#Fetch data: get the fasta files of the data and create data frames of the sequences

# Fetch data and create string sets of the fasta files 
recs_cytb <- entrez_fetch(db = "nucleotide", web_history = elephantidae_cytb$web_history, rettype = "fasta", retmax = elephantidae_cytb$count)
write(recs_cytb, "elephant_cytb.fasta")
cytb_stringSet <- readDNAStringSet("elephant_cytb.fasta")

recs_16S <- entrez_fetch(db = "nucleotide", web_history = elephantidae_16S$web_history, rettype = "fasta", retmax = elephantidae_16S$count)
write(recs_16S, "elephant_16S.fasta")
stringSet_16S <- readDNAStringSet("elephant_16S.fasta")

# Make data frames of the string sets
df_cytb <- data.frame(CytB_Title = names(cytb_stringSet), CytB_Sequence = paste(cytb_stringSet))
df_16S <- data.frame(Title_16S = names(stringSet_16S), Sequence_16S = paste(stringSet_16S))

# Add another column to the data frame that contains the species name of the data point. 
df_cytb$Species_Name <- word(df_cytb$CytB_Title, 2L, 3L)
df_cytb <- df_cytb[, c("CytB_Title", "Species_Name", "CytB_Sequence")]
View(df_cytb)

df_16S$Species_Name <- word(df_16S$Title_16S, 2L, 3L)
df_16S <- df_16S[, c("Title_16S", "Species_Name", "Sequence_16S")]
View(df_16S)

################################################################################
# Summarize the data

summary(df_16S)
summary(df_cytb)

# check to make sure there are no NAs in the sequence field for both data frames
df_16S %>%
  count(is.na(Sequence_16S))
df_cytb %>%
  count(is.na(CytB_Sequence))

# count the different instances of each species in both data frames
df_16S %>%
  group_by(Species_Name) %>%
  count(Species_Name)

df_cytb %>%
  group_by(Species_Name) %>%
  count(Species_Name)

# look at the summary for the length of sequences in each data frame
summary(nchar(df_16S$Sequence_16S))
summary(nchar(df_cytb$CytB_Sequence))

# create a histogram of the sequence lengths in each data frame
hist(nchar(df_cytb$CytB_Sequence), xlab = "Sequence Length", ylab = "
     Frequency", main = "Frequency Histogram of CytB Sequence Lengths")

hist(nchar(df_16S$Sequence_16S), xlab = "Sequence Length", ylab = "
     Frequency", main = "Frequency Histogram of 16S Sequence Lengths")

# Looking at the length of all the sequences in each set
df_16S %>%
  count(nchar(Sequence_16S))
df_16S_28 <- df_16S %>%
  filter(nchar(Sequence_16S) == 28)

write(recs_cytb, "elephant_cytb.fasta")
cytb_stringSet <- readDNAStringSet("elephant_cytb.fasta")


stringSet_16S

write(df_16S_28, "elephant_16S_28.fasta")

df_cytb %>%
  count(nchar(CytB_Sequence))

################################################################################
# Alignment
?muscle

cytb_alignment <- DNAStringSet(muscle::muscle(cytb_stringSet))
alignment_16S <- DNAStringSet(muscle::muscle(stringSet_16S))

BrowseSeqs(alignment_16S)
BrowseSeqs(cytb_alignment)

################################################################################

#############################################################################################

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