## 6210 Assignment 2

# The Differentiation of CytB sequences between Extinct and Extant Elephantidae species

# Load libraries 
library(rentrez)
library(BiocGenerics)
library(Biostrings)
library(tidyverse)
library(muscle)
library(DECIPHER)
library(ape)
library(randomForest)

#library(stringi)
#library(RSQLite)

elephant_search1 <- entrez_search(db = "nucleotide", term = "Elephantidae[ORGN] AND ND4[GENE] AND 3000:4000[SLEN]")



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
df_extinct$status <- "Extinct"
df_extant$status <- "Extant"

df_extinct1 <- df_extinct %>%
  mutate(Sequence0 = str_remove(Sequence, "^[-N]+")) %>%
  mutate(Sequence0 = str_remove(Sequence, "[-N]+$")) %>%
  mutate(Sequence0 = str_remove(Sequence, "-+")) %>%
  filter(str_count(Sequence0, "N") <= (0.05 * str_count(Sequence)))

df_extant1 <- df_extant %>%
  mutate(Sequence0 = str_remove(Sequence, "^[-N]+")) %>%
  mutate(Sequence0 = str_remove(Sequence, "[-N]+$")) %>%
  mutate(Sequence0 = str_remove(Sequence, "-+")) %>%
  filter(str_count(Sequence0, "N") <= (0.05 * str_count(Sequence)))

df_merge1 <- merge.data.frame(df_extant1, df_extinct1, all = T)
df_merge1$Sequence0 <- DNAStringSet(df_merge1$Sequence0)

df_merge1 <- cbind(df_merge1, as.data.frame(letterFrequency(df_merge1$Sequence0, letters = c("A", "C","G", "T"))))

df_merge1$Aprop <- (df_merge1$A) / (df_merge1$A + df_merge1$T + df_merge1$C + df_merge1$G)
df_merge1$Tprop <- (df_merge1$T) / (df_merge1$A + df_merge1$T + df_merge1$C + df_merge1$G)
df_merge1$Gprop <- (df_merge1$G) / (df_merge1$A + df_merge1$T + df_merge1$C + df_merge1$G)

df_merge1$Sequence0 <- as.character(df_merge1$Sequence0)

################################################################################

set.seed(16)
df_Validation <- df_merge1 %>%
  group_by(status) %>%
  sample_n(20)

set.seed(85)
df_Training <- df_merge1 %>%
  filter(!Title %in% df_Validation$Title) %>%
  group_by(status) %>%
  sample_n(60)

classifier <- randomForest::randomForest(x = df_Training[, 10:12], y = as.factor(df_Training$status), ntree = 50, importance = TRUE)

classifier

predict_Validation <- predict(classifier, df_Validation[, c(4, 10:12)])
predict_Validation

table(observed = df_Validation$status, predicted = predict_Validation)

################################################################################
df_extinct$status <- "Extinct"
df_extant$status <- "Extant"

merged <- merge(df_extant, df_extinct, all = T)
df_merged <- merged %>%
  mutate(Sequence_remove = str_remove_all(Sequence, "^N+|N+$|-")) %>%
  filter(str_count(Sequence_remove, "N") <= (0.01 * str_count(Sequence))) %>%
  filter(str_count(Sequence_remove) >= median(str_count(Sequence_remove)) - 75 & str_count(Sequence_remove) <= median(str_count(Sequence_remove)) + 75)


df_merged <- as.data.frame(df_merged)
df_merged$Sequence <- DNAStringSet(df_merged$Sequence)
df_merged$Sequence_remove <- DNAStringSet(df_merged$Sequence_remove)
names(df_merged$Sequence_remove) <- paste(df_merged$status, word(df_merged$Title, 1L), sep = " ")
df_merged$marker <- paste(df_merged$status, word(df_merged$Title, 1L), sep = " ")





merged_alignment <- DNAStringSet(muscle::muscle(df_merged$Sequence_remove))
BrowseSeqs(merged_alignment)

################################################################################
# Distance Matrix: calculate multiple distance matrices using different models
bin <- as.DNAbin(merged_alignment)

#Jukes and Canter
distanceMatrix_JC <- dist.dna(bin, model = "JC69", as.matrix = TRUE, pairwise.deletion = TRUE)
distanceMatrix_JC = as.dist(distanceMatrix_JC)
jc_cluster_single <- hclust(distanceMatrix_JC,method="single")
jc_cluster_complete <- hclust(distanceMatrix_JC,method="complete")

# Kimura's 2-parameter distance
distanceMatrix_Kim2P <- dist.dna(bin, model = "K80", as.matrix = TRUE, pairwise.deletion = TRUE)
distanceMatrix_Kim2P = as.dist(distanceMatrix_Kim2P)
Kim_cluster_single <- hclust(distanceMatrix_Kim2P,method="single")
Kim_cluster_complete <- hclust(distanceMatrix_Kim2P,method="complete")

distanceMatrix_Tamura <- dist.dna(bin, model = "T92", as.matrix = TRUE, pairwise.deletion = TRUE)
distanceMatrix_Tamura = as.dist(distanceMatrix_Tamura)

distanceMatrix_pariwise <- dist.dna(x = bin, model = "raw", as.matrix = T, pairwise.deletion = T)
distance_Matrix_pariwise <- as.dist(distanceMatrix_pariwise)
pairwise_cluster <- hclust(distance_Matrix_pariwise, method = "single")



plot(pairwise_cluster)
# Clustering


plot(as.phylo(jc_cluster_single), tip.color = c("red", "green", "blue")[cutree(jc_cluster_single, 3)])
plot(jc_cluster_complete)

################################################################################
df_merged$marker <- paste(df_merged$status, word(df_merged$Title, 1L), sep = " ")

df_merged_noPredict <- df_merged %>%
  filter(!Species_Name == "PREDICTED: Elephas", !Species_Name == "PREDICTED: Loxodonta", !marker== "Extinct FJ753551.1")

df_merged_noPredict <- as.data.frame(df_merged_noPredict)
df_merged_noPredict$Sequence <- DNAStringSet(df_merged_noPredict$Sequence)
df_merged_noPredict$Sequence_remove <- DNAStringSet(df_merged_noPredict$Sequence_remove)
names(df_merged_noPredict$Sequence_remove) <- paste(df_merged_noPredict$status, word(df_merged_noPredict$Title, 1L), sep = " ")

noPredict_alignment <- DNAStringSet(muscle::muscle(df_merged_noPredict$Sequence_remove))

bin_noPredict <- as.DNAbin(noPredict_alignment)

noPredict_distanceMatrix_pairwise <- dist.dna(x = bin_noPredict, model = "raw", as.matrix = T, pairwise.deletion = T)
View(noPredict_distanceMatrix_pairwise)


noPredict_distance_Matrix_pariwise <- as.dist(noPredict_distanceMatrix_pairwise)

noPredict_pairwise_cluster <- hclust(noPredict_distance_Matrix_pariwise, method = "single")


plot(as.phylo(noPredict_pairwise_cluster), tip.color = c("red", "blue")[cutree(noPredict_pairwise_cluster, 2)])

plot(as.phylo(noPredict_pairwise_cluster), tip.color = c("red", "blue", "green")[cutree(noPredict_pairwise_cluster, 3)])


cut <- cutree(noPredict_pairwise_cluster, 2)
View(cut)

cut

df_cut <- as.data.frame(cut)
df_cut %>%
  count(cut) 
cluster1 <- df_cut %>%
  filter(cut == 1)
cluster2 <- df_cut %>%
  filter(cut == 2)

df_cut$status <- word(row.names(df_cut), 1L, sep = " ")

df_cut %>% 
  group_by(cut) %>%
  count(status)

cut3 <- cutree(noPredict_pairwise_cluster, 3)

df_cut3 <- as.data.frame(cut3)

df_cut3$status <- word(row.names(df_cut3), 1L, sep = " ")

df_cut %>% 
  group_by(cut) %>%
  count(status)


################################################################################

cl <- stats::kmeans(distanceMatrix_JC, 2)

set.seed(96)
plot(distanceMatrix_JC, col = cl$cluster)








