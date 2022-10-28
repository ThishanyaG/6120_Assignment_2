###############################  6120 Assignment 2  ###############################

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


DataAcquisition <- function(search_terms, db="nucleotide", summary_term = "title") {
  search_results <- entrez_search(db=db, term = search_terms)
  search_results <- entrez_search(db=db, term = search_terms, retmax = search_results$count, use_history = T)

  summary <- entrez_summary(db = db, id = search_results$ids[1:300])
  View(extract_from_esummary(summary, summary_term))
  
  return(search_results)
}

Fasta_StringSet_DataFrame <- function(search_results, db = "nucleotide", file_name){

  recs <- entrez_fetch(db = db, web_history = search_results$web_history, rettype = "fasta", retmax = search_results$count)
  write(recs, file_name)
  string_set <- lapply(file_name, readDNAStringSet)
  
  # Make data frames of the string sets
  df <- data.frame(Title = names(string_set), Sequence = paste(string_set))
  
  # Add another column to the data frame that contains the species name of the data point. 
  df$Species_Name <- word(df$Title, 2L, 3L)
  df <- df[, c("Title", "Species_Name", "Sequence")]
  
  return(c(string_set, df))
}

ClassifierPrep <- function(df){
  
  df %>%
  mutate(Sequence0 = str_remove(Sequence, "^[-N]+")) %>%
  mutate(Sequence0 = str_remove(Sequence, "[-N]+$")) %>%
  mutate(Sequence0 = str_remove(Sequence, "-+")) %>%
  filter(str_count(Sequence0, "N") <= (0.05 * str_count(Sequence))) 
  
  df$Sequence0 <- DNAStringSet(df$Sequence0)
  df <- cbind(df, as.data.frame(letterFrequency(df$Sequence0, letters = c("A", "C","G", "T"))))
  
  df$Aprop <- (df$A) / (df$A + df$T + df$C + df$G)
  df$Tprop <- (df$T) / (df$A + df$T + df$C + df$G)
  df$Gprop <- (df$G) / (df$A + df$T + df$C + df$G)
  
  df$Sequence0 <- as.character(df$Sequence0)
  
  return(df)
}

RandomForestClassifier <- function(df, group, seed1=16, sample1=20, seed2=85, sample2=60, n_tree=50){
  
  set.seed(seed1)
  df_Validation <- df %>%
    group_by(group) %>%
    sample_n(sample1)
  
  set.seed(seed2)
  df_Training <- df %>%
    filter(!Title %in% df_Validation$Title) %>%
    group_by(group) %>%
    sample_n(sample2)
  
  classifier <- randomForest::randomForest(x = df_Training[, 10:12], y = as.factor(df_Training$group), ntree = ntree, importance = TRUE)
  
  predict_Validation <- predict(classifier, df_Validation[, c(4, 10:12)])
  
  return(classifier, predict_Validation)
}

ClusteringPrep <- function(df){
  df %>%
    mutate(Sequence_remove = str_remove_all(Sequence, "^N+|N+$|-")) %>%
    filter(str_count(Sequence_remove, "N") <= (0.01 * str_count(Sequence))) %>%
    filter(str_count(Sequence_remove) >= median(str_count(Sequence_remove)) - 75 & str_count(Sequence_remove) <= median(str_count(Sequence_remove)) + 75)
  
  df <- as.data.frame(df)
  df$Sequence <- DNAStringSet(df$Sequence)
  df$Sequence_remove <- DNAStringSet(df$Sequence_remove)
  names(df$Sequence_remove) <- paste(df$status, word(df$Title, 1L), sep = " ")
  df$marker <- paste(df$status, word(df$Title, 1L), sep = " ")
  
  return(df)
}

DistanceMatrix1 <- function(alignment, model, cluster_method){
  bin <- as.DNAbin(alignment)
  
  distance_matrix <- dist.dna(bin, model = mpdel, as.matrix = TRUE, pairwise.deletion = TRUE)
  distance_matrix = as.dist(distance_matrix)
  cluster <- hclust(distance_matrix, method = cluster_method)
  
}



extinct_search_terms = "Elephas antiquus[ORGN] OR Mammut[ORGN] OR Mammuthus[ORGN] OR Elephas cypriotes[ORGN] OR Elephas maximus asurus[ORGN] AND cytochrome b AND 400:1000[SLEN]"
extinct_elephants <- DataAcquisition(extinct_search_terms)
c(extinct_string_set, df_extinct) <- Fasta_StringSet_DataFrame(extinct_elephants, file_name = "extinct_elephant_cytb.fasta")
  
  
  
extant_search_terms = "Loxodonta africana[ORGN] OR Loxodonta cyclotis OR (Elephas maximus[ORGN] NOT Elephas maximus asurus[ORGN]) AND cytochrome b AND 500:1000[SLEN]"
extant_elephants <- DataAcquisition(extant_search_terms)
  
?apply
