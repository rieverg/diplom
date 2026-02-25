.libPaths(c("/home/mdashian/R/x86_64-redhat-linux-gnu-library/4.0", .libPaths()))

library(Biostrings)
library(BiocGenerics)
library(GenomeInfoDb)
library(ape)
library(dplyr)
library(stringr)
library(msa)



setwd("/home/mdashian/diplom/ebola")
input_directory <- "/home/mdashian/diplom/ebola"
fragment_lengths <- c(100, 500, 1000)
output_dir <- "/home/mdashian/diplom/ebola/out"
n_total <- 100

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

fasta_files <- list.files(path = input_directory, pattern = "\\.fasta$|\\.fa$|\\.fna$", full.names = TRUE)

all_sequences <- DNAStringSet()
for (file in fasta_files) {
  seqs <- readDNAStringSet(file)
  min_length <- max(fragment_lengths)
  seqs <- seqs[width(seqs) >= min_length] 
  all_sequences <- c(all_sequences, seqs)
}

valid_seqs <- DNAStringSet(use.names=TRUE)
#for (i in 1:length(all_sequences)) {
#  seq <- all_sequences[i]
#  seq_str <- as.character(seq)
  
#  if (!all(strsplit(seq_str, "")[[1]] %in% c("a", "t", "c", "g", "-"))) next
  
#  seq <- DNAStringSet(seq)
#  valid_seqs <- append(valid_seqs, seq)
#}

#filtered <- DNAStringSet(valid_seqs, use.names=TRUE) не знаю что сделать с уже выровненными, мб сначала контроль потом выравнивание

filtered = all_sequences


options(error = function() {
  traceback(3)
  quit(status = 1)
})

message("Starting at: ", Sys.time())
message("Filtered sequences: ", length(filtered))


fast_similarity <- function(seq1, seq2) {
  s1 <- as.character(seq1)
  s2 <- as.character(seq2)
  
  if (nchar(s1) > nchar(s2)) {
    temp <- s1
    s1 <- s2
    s2 <- temp
  }
  
  eq <- 0  # совпадения
  neq <- 0 # несовпадения
  
  for (i in 1:nchar(s1)) {
    char1 <- substr(s1, i, i)
    char2 <- substr(s2, i, i)
    
    if (char1 == "-" || char2 == "-" || char1 == "N" || char2 == "N") next
    
    if (char1 == char2) {
      eq <- eq + 1
    } else {
      neq <- neq + 1
    }
  }
  
  if (neq + eq == 0) return(0)
  return(100 * eq / (neq + eq))  
}


# 1. Создание выборки "случайные + родственники"
set.seed(123)
n_seeds <- 5
identity_threshold <- 95
seeds <- sample(1:length(filtered), n_seeds)

relatives <- list()
for (seed in seeds) {
  seq <- filtered[seed]
  pids <- sapply(1:length(filtered), function(i) {
    if (i == seed) return(100)
    aln <- pairwiseAlignment(seq, filtered[i], type="global")
    pid(aln, type="PID1")
  })
  relatives[[as.character(seed)]] <- which(pids >= identity_threshold)
}

all_relatives <- unique(unlist(relatives))

if (length(all_relatives) > n_total) {
  selected_random <- sample(all_relatives, n_total)
} else {
  selected_random <- all_relatives
}

random_set <- filtered[selected_random]
writeXStringSet(random_set, file.path(output_dir, "random_plus_relatives.fasta"))


pairs_to_check <- min(1000, length(filtered))
max_similarity <- -Inf
best_pair <- c(1, 2)

for (i in 1:(pairs_to_check-1)) {
  for (j in (i+1):pairs_to_check) {
    current_sim <- fast_similarity(filtered[i], filtered[j])
    if (current_sim > max_similarity) {
      max_similarity <- current_sim
      best_pair <- c(i, j)
    }
  }
}

selected_similar <- best_pair


selected_diverse <- sample(1:length(filtered), 2)

while (length(selected_diverse) < n_total) {
  min_pids <- sapply(setdiff(1:length(filtered), selected_diverse), function(i) {
    min(sapply(selected_diverse, function(j) {
      fast_similarity(filtered[i], filtered[j])
    }))
  })
  
  next_seq <- setdiff(1:length(filtered), selected_diverse)[which.min(min_pids)]
  selected_diverse <- c(selected_diverse, next_seq)
  
  
  }


diverse_set <- filtered[selected_diverse]
writeXStringSet(diverse_set, file.path(output_dir, "most_diverse.fasta"))


message("Все выборки созданы и сохранены в: ", output_dir)
