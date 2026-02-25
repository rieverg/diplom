.libPaths(c("/home/mdashian/R/x86_64-redhat-linux-gnu-library/4.0", .libPaths()))

library(Biostrings)
library(BiocGenerics)
library(GenomeInfoDb)
library(ape)
library(dplyr)
library(stringr)
library(msa)



# Установка рабочей директории и параметров
setwd("/home/mdashian/diplom/morbillivirus_cluster")
input_directory <- "/home/mdashian/diplom/morbillivirus_cluster"
fragment_lengths <- c(100, 500, 1000)
output_dir <- "/home/mdashian/diplom/morbillivirus_cluster"
n_total <- 100

# Создаем выходную директорию
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Загрузка и фильтрация последовательностей
fasta_files <- "/home/mdashian/diplom/morbillivirus_cluster/1_aligned.fasta"

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
#filtered = filtered[(1:30)]
filtered = all_sequences



# Включите подробное логирование
options(error = function() {
  traceback(3)
  quit(status = 1)
})

message("Starting at: ", Sys.time())
message("Filtered sequences: ", length(filtered))


# Быстрая функция сравнения - возвращает процент СХОДСТВА (как pid)
fast_similarity <- function(seq1, seq2) {
  # Преобразуем в character для быстрого доступа
  s1 <- as.character(seq1)
  s2 <- as.character(seq2)
  
  # Делаем s1 более короткой последовательностью
  if (nchar(s1) > nchar(s2)) {
    temp <- s1
    s1 <- s2
    s2 <- temp
  }
  
  eq <- 0  # совпадения
  neq <- 0 # несовпадения
  
  # Быстрое посимвольное сравнение
  for (i in 1:nchar(s1)) {
    char1 <- substr(s1, i, i)
    char2 <- substr(s2, i, i)
    
    # Пропускаем gaps и N
    if (char1 == "-" || char2 == "-" || char1 == "N" || char2 == "N") next
    
    if (char1 == char2) {
      eq <- eq + 1
    } else {
      neq <- neq + 1
    }
  }
  
  # Возвращаем процент СХОДСТВА (как pid), а не различий!
  if (neq + eq == 0) return(0)
  return(100 * eq / (neq + eq))  # процент идентичности
}




# 3. Создание выборки "наиболее похожие"
# Используем вашу быструю функцию сравнения

# Поиск двух самых похожих последовательностей
pairs_to_check <- min(1000, length(filtered))
max_similarity <- -Inf
best_pair <- c(1, 2)

for (i in 1:(pairs_to_check-1)) {
  for (j in (i+1):pairs_to_check) {
    # ЗАМЕНА: fast_similarity вместо pid(pairwiseAlignment())
    current_sim <- fast_similarity(filtered[i], filtered[j])
    if (current_sim > max_similarity) {
      max_similarity <- current_sim
      best_pair <- c(i, j)
    }
  }
}

selected_similar <- best_pair

# Добавление остальных похожих последовательностей
while (length(selected_similar) < n_total) {
  avg_similarities <- sapply(setdiff(1:length(filtered), selected_similar), function(i) {
    mean(sapply(selected_similar, function(j) {
      # ЗАМЕНА: fast_similarity вместо pid(pairwiseAlignment())
      fast_similarity(filtered[i], filtered[j])
    }))
  })
  
  # Выбираем последовательность с максимальным средним сходством
  next_seq <- setdiff(1:length(filtered), selected_similar)[which.max(avg_similarities)]
  selected_similar <- c(selected_similar, next_seq)
  
  
  }


similar_set <- filtered[selected_similar]
writeXStringSet(similar_set, file.path(output_dir, "most_similar_1.fasta"))


message("Все выборки созданы и сохранены в: ", output_dir)