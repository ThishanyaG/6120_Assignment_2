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
library(phangorn)
library(cluster)
library(vegan)


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
df_extinct$status <- "Extinct"

View(df_extinct)

# And repeat for extant species
recs_extant <- entrez_fetch(db = "nucleotide", web_history = extant_elephants$web_history, rettype = "fasta", retmax = extant_elephants$count)
write(recs_extant, "extant_elephant_cytb.fasta")
extant_cytb_stringSet <- readDNAStringSet("extant_elephant_cytb.fasta")

df_extant <- data.frame(Title = names(extant_cytb_stringSet), Sequence = paste(extant_cytb_stringSet))

# Add another column to the data frame that contains the species name of the data point. 
df_extant$Species_Name <- word(df_extant$Title, 2L, 3L)
df_extant <- df_extant[, c("Title", "Species_Name", "Sequence")]
df_extant$status <- "Extant"
View(df_extant)

rm(recs_extant, recs_extinct)
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
par(mfrow=c(2,1))

extinct_length <- hist(nchar(df_extinct$Sequence), xlab = "Sequence Length", ylab = "Frequency", main = "Frequency Histogram of CytB Sequence Lengths in Extinct Species")

extant_length <- hist(nchar(df_extant$Sequence), xlab = "Sequence Length", ylab = " Frequency", main = "Frequency Histogram of CytB Sequence Lengths in Extant Species")


# Looking at the length of all the sequences in each set
df_extinct %>%
  count(nchar(Sequence))

df_extant %>%
  count(nchar(Sequence))

# Alignment: do initial alignments for both data frames, separately for now. 
extinct_alignment <- DNAStringSet(muscle::muscle(extinct_cytb_stringSet, maxiters = 2))
extant_alignment <- DNAStringSet(muscle::muscle(extant_cytb_stringSet, maxiters = 2))

BrowseSeqs(extinct_alignment)
BrowseSeqs(extant_alignment)

rm(extant_alignment, extinct_alignment)
################################################################################
# Prepare data for clustering

# merge the 2 data frames into 1. remove any leading or trailing N's and gaps from the data, then filter the  data frame to remove all sequences that have a high N counts (more that 1%) or sequences with large size variation. The edited sequences are added to a new column.
merged <- merge(df_extant, df_extinct, all = T)
df_merged <- merged %>%
  mutate(Sequence_remove = str_remove_all(Sequence, "^N+|N+$|-")) %>%
  filter(str_count(Sequence_remove, "N") <= (0.01 * str_count(Sequence))) %>%
  filter(str_count(Sequence_remove) >= median(str_count(Sequence_remove)) - 75 & str_count(Sequence_remove) <= median(str_count(Sequence_remove)) + 75)

# Ensure the new data frame is set as a data frame and make the 2 sequence columns into DNA String set objects. The names of the edited sequences are set to the status and accession number to ensure unique identifiers. A new columns is also added that contains this identifier
df_merged <- as.data.frame(df_merged)
df_merged$Sequence <- DNAStringSet(df_merged$Sequence)
df_merged$Sequence_remove <- DNAStringSet(df_merged$Sequence_remove)
names(df_merged$Sequence_remove) <- paste(df_merged$status, word(df_merged$Title, 1L), sep = " ")
df_merged$marker <- paste(df_merged$status, word(df_merged$Title, 1L), sep = " ")

# do an initial alignment of the newly merged data
merged_alignment <- DNAStringSet(muscle::muscle(df_merged$Sequence_remove, maxiters = 2, log = "log.txt"))
BrowseSeqs(merged_alignment)

# do an alignment with a higher gap penalty
merged_alignment_gap <- DNAStringSet(muscle::muscle(df_merged$Sequence_remove, maxiters = 2, gapopen = -1000))
BrowseSeqs(merged_alignment_gap)

rm(merged_alignment, merged_alignment_gap)

################################################################################
# After analyzing the alignments, it is clear that there are a few sequences that vary wildly from the other sequences. Thus, they will be removed and added to a new data frame for further analysis
df_merged$Sequence <- as.character(df_merged$Sequence)
df_merged$Sequence_remove <- as.character(df_merged$Sequence_remove)

df_merged_noPredict <- df_merged %>%
  filter(!Species_Name == "PREDICTED: Elephas", !Species_Name == "PREDICTED: Loxodonta", !marker== "Extinct FJ753551.1")

# Ensure the new data frame is set as a data frame and make the 2 sequence columns into DNA String set objects. The names of the edited sequences are set to the status and accession number to ensure unique identifiers. A new columns is also added that contains this identifier
df_merged_noPredict$Sequence <- DNAStringSet(df_merged_noPredict$Sequence)
df_merged_noPredict$Sequence_remove <- DNAStringSet(df_merged_noPredict$Sequence_remove)
df_merged_noPredict <- as.data.frame(df_merged_noPredict)
names(df_merged_noPredict$Sequence_remove) <- paste(df_merged_noPredict$status, word(df_merged_noPredict$Title, 1L), sep = " ")

# The merged data is aligned
noPredict_alignment <- DNAStringSet(muscle::muscle(df_merged_noPredict$Sequence_remove))
BrowseSeqs(noPredict_alignment)

################################################################################
#Distance Matrix:

# create a DNA bin object from the previous alignment
bin_noPredict <- as.DNAbin(noPredict_alignment)

# Calculate the distance matrix with pairwise deletion
noPredict_distanceMatrix_pairwise <- dist.dna(x = bin_noPredict, model = "raw", as.matrix = T, pairwise.deletion = T)
noPredict_distance_Matrix_pairwise <- as.dist(noPredict_distanceMatrix_pairwise)

################################################################################
# Clustering: 

# Complete clustering method:
noPredict_pairwise_cluster <- hclust(noPredict_distance_Matrix_pairwise, method = "complete")

# Plot the results of the clustering
plot(as.phylo(noPredict_pairwise_cluster), tip.color = c("red", "blue", "green")[cutree(noPredict_pairwise_cluster, 3)])

plot(as.phylo(noPredict_pairwise_cluster), show.tip.label = FALSE, main = "Complete Clustering Tree", tip)

plot(as.dendrogram(noPredict_pairwise_cluster), nodePar = list(col = c("red", "blue", "green")[cutree(noPredict_pairwise_cluster, 3)]), leaflab = "none")
??nodePar
# Analyse the tips of the trees based on how they have been clustered. We cut the tree based off of the outermost branches and summerise
cut <- cutree(noPredict_pairwise_cluster, 2)
df_cut <- as.data.frame(cut)
df_cut$status <- word(row.names(df_cut), 1L, sep = " ")
df_cut %>% 
  group_by(cut) %>%
  count(status)

cut3 <- cutree(noPredict_pairwise_cluster, 3)
df_cut3 <- as.data.frame(cut3)
df_cut3$status <- word(row.names(df_cut3), 1L, sep = " ")
df_cut3 %>% 
  group_by(cut3) %>%
  count(status)

cut4 <- cutree(noPredict_pairwise_cluster, 4)
df_cut4 <- as.data.frame(cut4)
df_cut4$status <- word(row.names(df_cut4), 1L, sep = " ")
df_cut4 %>% 
  group_by(cut4) %>%
  count(status)

################################################################################
#K-Means:
# Use a function from Nallathambi, J. (2018) from https://medium.com/codesmart/r-series-k-means-clustering-silhouette-794774b46586 to Calculate multiple silhouette scores and plot the scores to find the ideal number of clusters to perform K-means clustering. According to the plot, k=3 is ideal
silhouette_score <- function(k){
  km <- kmeans(noPredict_distanceMatrix_pairwise, centers = k)
  ss <- silhouette(km$cluster, dist(noPredict_distanceMatrix_pairwise))
  mean(ss[,3])
}

k <- 2:10
sil_func <- lapply(k, silhouette_score)
plot(k, type='b', sil_func, main ="Silhouette Scores at Mulitple K-values",xlab = "Value of K", ylab = "Silhouette Score")

# Perform the K-means with k=3, calculate the silhouette score and plot the silhouette index
cl3 <- kmeans(noPredict_distance_Matrix_pairwise, 3)
sil_clus3 <- silhouette(cl3$cluster, dist(noPredict_distance_Matrix_pairwise))

grDevices::windows()
plot(sil_clus3, main = "3 Cluster Silhouette Plot")
plot(noPredict_distance_Matrix_pairwise, col = cl3$cluster)

