.libPaths(c("/home/mdashian/R/x86_64-redhat-linux-gnu-library/4.0", .libPaths()))

required_packages <- c("later", "TreeDist", "Biostrings", "msa", "ape", "phangorn", "dplyr", "stringr", "foreach", "ggplot2")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE, repos = "https://cloud.r-project.org")
  }
  library(pkg, character.only = TRUE)
}
library(Biostrings)
library(msa)
library(ape)
library(phangorn)
library(dplyr)
library(stringr)
library(TreeDist)
library(ggplot2)
library(parallel)

input_directory <- "/home/mdashian/diplom/syncytial_cluster/"
sample_size <- 100




#fasta_files <- "/home/mdashian/diplom/sequences/cluster/most_similar_4_aligned.fasta"
fasta_files <- "/home/mdashian/diplom/syncytial_cluster/most_similar_aligned.fasta"


all_sequences <- DNAStringSet()
for (file in fasta_files) {
  seqs <- readDNAStringSet(file)
  min_length <- 10000
  seqs <- seqs[width(seqs) >= min_length] 
  if (length(seqs) > sample_size) {
    seqs <- seqs[sample(1:length(seqs), sample_size)]
  }
  all_sequences <- c(all_sequences, seqs)
}
max_length <- max(width(all_sequences))

filtered=all_sequences




new=as.DNAbin(filtered)

full_aligned <- as.phyDat(new)

dm <- dist.ml(full_aligned)
full_tree <- NJ(dm)





fit <- pml(full_tree, data = full_aligned, model = "GTR") #обязательно проверить для нового вируса лучшую модель!!!!!
fit_optimized <- optim.pml(fit,
                           model = "GTR",
                           optInv = TRUE,
                           optGamma = TRUE,
                           optNni = TRUE,
                           rearrangement = "NNI", control = pml.control(epsilon = 1e-10)) 

bs_trees_full <- bootstrap.pml(fit_optimized, bs = 100, optNni = TRUE,
                         multicore = TRUE, mc.cores = detectCores()-1) #выключить или включить nni в обоих функциях!!!


# Добавление поддержки на дерево
fit_optimized$tree <- plotBS(fit_optimized$tree, bs_trees_full, type = "none")



seq_length <- width(seqs[1])


fragment_lengths <- c(500, 1000, 2000, 3000, 3500, 4000, 4500, 5000, 6000, 8000, 10000)
rf_distances <- numeric(length(fragment_lengths))


build_fragment_tree <- function(seqs, fragment_len, full_tree) {
  fragments <- DNAStringSet()
  

  for (i in 1:length(seqs)) {
    seq <- seqs[i]
    seq_length <- width(seq)
    start <- 2500
    if (seq_length-start > fragment_len) {
      fragment <- subseq(seq, start, start + fragment_len - 1)
      fragments <- c(fragments, DNAStringSet(fragment))
    } else {
    fragment <- subseq(seq, seq_length-start+1, seq_length)
    fragments <- c(fragments, DNAStringSet(fragment))
    }
  }
  
    
 
  fragments_bin <- as.DNAbin(fragments)
  aligned_frag_phy <- as.phyDat(fragments_bin)
  

  dm_frag <- dist.ml(aligned_frag_phy)
  frag_tree <- NJ(dm_frag)
  

  fit_frag <- pml(frag_tree, data = aligned_frag_phy, model = "GTR")
  fit_optimized_frag <- optim.pml(fit_frag,
                                 model = "GTR",
                                 optInv = TRUE,
                                 optGamma = TRUE,
                                 optNni = TRUE,
                                 rearrangement = "NNI",  control = pml.control(epsilon = 1e-10))
  
  bs_trees_frag <- bootstrap.pml(fit_optimized_frag, bs = 100, optNni = TRUE,
                         multicore = TRUE, mc.cores = detectCores()-1)
  


 fit_optimized_frag$tree <- plotBS(fit_optimized_frag$tree, bs_trees_frag, type = "none")
  
  
  rf_dist <- RF.dist(fit_optimized_frag$tree, fit_optimized$tree, 
                    normalize = TRUE, check.labels = TRUE, rooted = FALSE)
  gen_dist = MutualClusteringInfo(
  fit_optimized$tree,
  fit_optimized_frag$tree,
  normalize = TRUE,
  reportMatching = FALSE,
  diag = FALSE
)
  result=data.frame(gen_dist, rf_dist)
  return(result)
  
  
}


all_results <- data.frame()

for (j in seq_along(fragment_lengths)) {
  cat("Анализируем длину фрагмента:", fragment_lengths[j], "bp\n")
  
  n_replicates <- 1
  replicate_results <- data.frame()
  
  for (rep in 1:n_replicates) {
    cat("  Повторение", rep, "из", n_replicates, "\n")
    result <- build_fragment_tree(filtered, fragment_lengths[j], fit_optimized$tree)
    
   
    result_row <- data.frame(
      Fragment_Length = fragment_lengths[j],
      Replicate = rep,
      RF_Distance = result$rf_dist,
      Gen_Distance = result$gen_dist
    )
    
    replicate_results <- rbind(replicate_results, result_row)
  }
  
  all_results <- rbind(all_results, replicate_results)
  cat("  Среднее RF расстояние:", mean(replicate_results$RF_Distance), "\n")
  cat("  Среднее Gen расстояние:", mean(replicate_results$Gen_Distance), "\n\n")
  
  
  
}



results_summary <- all_results %>%
  group_by(Fragment_Length) %>%
  summarise(
    RF_Distance_Mean = mean(RF_Distance, na.rm = TRUE),
    RF_Distance_SD = sd(RF_Distance, na.rm = TRUE),
    Gen_Distance_Mean = mean(1-Gen_Distance, na.rm = TRUE),
    Gen_Distance_SD = sd(1-Gen_Distance, na.rm = TRUE)
  ) %>%
  ungroup()
  
print(all_results)
print(results_summary)
  
  

ggplot(results_summary, aes(x = Fragment_Length, y = Gen_Distance_Mean)) +
  geom_point(size = 3, color = "darkred") +
  geom_line(color = "darkred", linewidth = 1) +
  geom_errorbar(aes(ymin = Gen_Distance_Mean - Gen_Distance_SD, 
                    ymax = Gen_Distance_Mean + Gen_Distance_SD), 
                width = 50, color = "darkred", alpha = 0.6) +
  labs(title = "Зависимость Generalised RF расстояния от длины фрагмента",
       x = "Длина фрагмента (bp)",
       y = "Generalised расстояние") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("gen_distance_plot_similar_2500.png", width = 8, height = 6, dpi = 300)



library(tidyr)
results_long <- results_summary %>%
  select(Fragment_Length, RF_Distance_Mean, Gen_Distance_Mean) %>%
  pivot_longer(cols = c(RF_Distance_Mean, Gen_Distance_Mean), 
               names_to = "Metric", 
               values_to = "Distance")


ggplot(results_long, aes(x = Fragment_Length, y = Distance, color = Metric, group = Metric)) +
  geom_point(size = 3) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = c("RF_Distance_Mean" = "steelblue", 
                               "Gen_Distance_Mean" = "darkred"),
                     labels = c("gen RF Distance", "RF Distance")) +
  labs(title = "Сравнение RF и Gen RF расстояний",
       x = "Длина фрагмента (bp)",
       y = "Расстояние",
       color = "Метрика") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "bottom")

ggsave("combined_distance_plot_similar_2500.png", width = 8, height = 6, dpi = 300)
