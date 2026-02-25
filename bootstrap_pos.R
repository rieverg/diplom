.libPaths(c("/home/mdashian/R/x86_64-redhat-linux-gnu-library/4.0", .libPaths()))

library(Biostrings)
library(msa)
library(ape)
library(phangorn)
library(dplyr)
library(stringr)
library(foreach)
library(parallel)
library(TreeDist)

# Создаем директорию для сохранения деревьев
tree_dir <- "saved_trees_morb"
sample_size <- 100

fasta_files <- "/home/mdashian/diplom/morbillivirus_cluster/most_similar_aligned.fasta"

all_sequences <- DNAStringSet()
for (file in fasta_files) {
  seqs <- readDNAStringSet(file)
  min_length <- 15000
  seqs <- seqs[width(seqs) >= min_length] 
  if (length(seqs) > sample_size) {
    seqs <- seqs[sample(1:length(seqs), sample_size)]
  }
  all_sequences <- c(all_sequences, seqs)
}

filtered <- all_sequences

# Построение полногеномного дерева
new <- as.DNAbin(filtered)
full_aligned <- as.phyDat(new)

fit <- pml(NJ(dist.ml(full_aligned)), data = full_aligned, model = "GTR")

fit_optimized <- optim.pml(fit,
                           model = "GTR",
                           optInv = TRUE,
                           optGamma = TRUE,
                           optNni = TRUE,
                           rearrangement = "NNI", 
                           control = pml.control(epsilon = 1e-10))

# Бутстреп для полного дерева
bs_trees_full <- bootstrap.pml(fit_optimized, bs = 100, optNni = TRUE,
                               multicore = TRUE, mc.cores = detectCores()-1)
fit_optimized$tree <- plotBS(fit_optimized$tree, bs_trees_full, type = "none")

# Задаем диапазон длин фрагментов для анализа
fragment_lengths <- c(3000)

# Функция для построения дерева для фрагмента с заданным началом и длиной
build_fragment_tree <- function(seqs, fragment_len, start_pos, full_tree) {
  write.tree(fit_optimized$tree, file.path(tree_dir, "full_genome_tree.tre"))
  fragments <- DNAStringSet()
  
  # Выбираем фрагменты с заданным началом и длиной
  for (i in 1:length(seqs)) {
    seq <- seqs[i]
    seq_length <- width(seq)
    
    # Если фрагмент не помещается от текущей стартовой позиции, берем его от конца
    if (start_pos + fragment_len - 1 <= seq_length) {
      # Фрагмент помещается нормально
      fragment <- subseq(seq, start_pos, start_pos + fragment_len - 1)
    } else {
      # Фрагмент не помещается - берем от конца последовательности
      end_pos <- seq_length
      start_pos_adj <- end_pos - fragment_len + 1
      # Исправлено: добавлено извлечение фрагмента
      fragment <- subseq(seq, start_pos_adj, end_pos)
    }
    fragments <- c(fragments, DNAStringSet(fragment))
  }
  
  # Проверяем, что осталось достаточно последовательностей
  if (length(fragments) < 3) {
    warning("Недостаточно последовательностей для построения дерева")
    return(NA)
  }
  
  # Убедимся, что все фрагменты имеют одинаковую длину
  if (length(unique(width(fragments))) > 1) {
    max_frag_len <- max(width(fragments))
    fragments <- DNAStringSet(sapply(fragments, function(x) {
      if (width(x) < max_frag_len) {
        return(DNAStringSet(paste0(as.character(x), 
                                 paste(rep("-", max_frag_len - width(x)), collapse = ""))))
      } else {
        return(x)
      }
    }))
  }
  
  # Преобразуем в формат для филогенетического анализа
  fragments_bin <- as.DNAbin(fragments)
  aligned_frag_phy <- as.phyDat(fragments_bin)
  
  # Построение дерева для фрагментов
  dm_frag <- dist.ml(aligned_frag_phy)
  frag_tree <- NJ(dm_frag)
  
  # Оптимизация дерева
  fit_frag <- pml(frag_tree, data = aligned_frag_phy, model = "GTR")
  fit_optimized_frag <- optim.pml(fit_frag,
                                  model = "GTR",
                                  optInv = TRUE,
                                  optGamma = TRUE,
                                  optNni = TRUE,
                                  rearrangement = "NNI",
                                  control = pml.control(epsilon = 1e-10))
  
  # Бутстреп для фрагментного дерева
  bs_trees_frag <- bootstrap.pml(fit_optimized_frag, bs = 100, optNni = TRUE,
                                 multicore = TRUE, mc.cores = detectCores()-1)

  fit_optimized_frag$tree <- plotBS(fit_optimized_frag$tree, bs_trees_frag, type = "none")
  
  # Вычисляем RF расстояние
  rf_dist <- RF.dist(fit_optimized_frag$tree, full_tree, 
                     normalize = TRUE, check.labels = TRUE, rooted = FALSE)
  
  # Вычисляем Mutual Clustering Info
  gen_dist <- MutualClusteringInfo(
    full_tree,
    fit_optimized_frag$tree,
    normalize = TRUE,
    reportMatching = FALSE,
    diag = FALSE
  )
  
  # Возвращаем оба показателя
  result <- data.frame(rf_distance = rf_dist, mci = gen_dist)
  return(result)
}  # ЗАКРЫВАЕМ ФУНКЦИЮ build_fragment_tree ЗДЕСЬ

# Функция для анализа всех возможных начальных позиций
analyze_all_start_positions <- function(seqs, fragment_len, full_tree, step_size = 1000) {
  # Определяем минимальную длину последовательности
  min_seq_length <- min(width(seqs))
  max_seq_length <- max(width(seqs))
  
  cat("  Минимальная длина последовательности:", min_seq_length, "\n")
  cat("  Максимальная длина последовательности:", max_seq_length, "\n")
  cat("  Длина фрагмента:", fragment_len, "\n")
  
  # Определяем возможные начальные позиции
  max_start_pos <- 16000
  
  # Создаем последовательность начальных позиций с заданным шагом
  start_positions <- seq(1, max_start_pos, by = step_size)
  
  # Добавляем последнюю возможную позицию
  if (tail(start_positions, 1) != max_start_pos) {
    start_positions <- c(start_positions, max_start_pos)
  }
  
  results <- data.frame(start_position = integer(), 
                       rf_distance = numeric(),
                       mci = numeric())
  
  # Перебираем все начальные позиции
  for (start_pos in start_positions) {
    cat("  Анализируем начальную позицию:", start_pos, "\n")
    
    tree_result <- build_fragment_tree(seqs, fragment_len, start_pos, full_tree)
    
    if (!is.na(tree_result$rf_distance) && !is.na(tree_result$mci)) {
      results <- rbind(results, data.frame(
        start_position = start_pos, 
        rf_distance = tree_result$rf_distance,
        mci = tree_result$mci
      ))
    }
  }
  
  return(results)
}

# Основной цикл анализа для каждой длины фрагмента
all_results <- data.frame()

for (j in seq_along(fragment_lengths)) {
  cat("Анализируем длину фрагмента:", fragment_lengths[j], "bp\n")
  
  # Анализируем все возможные начальные позиции
  fragment_results <- analyze_all_start_positions(filtered, fragment_lengths[j], fit_optimized$tree, step_size = 1000)
  
  if (nrow(fragment_results) > 0) {
    fragment_results$fragment_length <- fragment_lengths[j]
    all_results <- rbind(all_results, fragment_results)
    
    # Выводим статистику по текущей длине фрагмента
    cat("  RF расстояние:\n")
    cat("    Среднее:", mean(fragment_results$rf_distance, na.rm = TRUE), "\n")
    cat("    Минимальное:", min(fragment_results$rf_distance, na.rm = TRUE), "\n")
    cat("    Максимальное:", max(fragment_results$rf_distance, na.rm = TRUE), "\n")
    
    cat("  Mutual Clustering Info:\n")
    cat("    Среднее:", mean(fragment_results$mci, na.rm = TRUE), "\n")
    cat("    Минимальное:", min(fragment_results$mci, na.rm = TRUE), "\n")
    cat("    Максимальное:", max(fragment_results$mci, na.rm = TRUE), "\n")
    
    cat("  Количество проанализированных позиций:", nrow(fragment_results), "\n\n")
  }
}

# Визуализация результатов
library(ggplot2)

if (nrow(all_results) > 0) {
  # График 1: RF расстояние vs начальная позиция
  p1 <- ggplot(all_results, aes(x = start_position, y = rf_distance)) +
    geom_point(size = 2, color = "steelblue") +
    geom_line(color = "steelblue", linewidth = 1) +
    labs(title = paste("Зависимость RF расстояния от начальной позиции фрагмента"),
         subtitle = paste("Длина фрагмента:", unique(all_results$fragment_length), "bp"),
         x = "Начальная позиция фрагмента",
         y = "RF расстояние (нормированное)") +
    theme_minimal()
  
  print(p1)
  ggsave("rf_distance_vs_start_position_similar_morb.png", width = 10, height = 6, dpi = 300)
  
  # График 2: MCI vs начальная позиция
  p2 <- ggplot(all_results, aes(x = start_position, y = (1-mci))) +
    geom_point(size = 2, color = "darkgreen") +
    geom_line(color = "darkgreen", linewidth = 1) +
    geom_smooth() +
    labs(title = paste("Зависимость Generalised RF distance от начальной позиции фрагмента"),
         subtitle = paste("Длина фрагмента:", unique(all_results$fragment_length), "bp"),
         x = "Начальная позиция фрагмента",
         y = "Generalised RF distance (нормированное)") +
    theme_minimal()
  
  print(p2)
  ggsave("mci_vs_start_position_similar_morb.png", width = 10, height = 6, dpi = 300)
  
  # График 3: Сравнение RF и MCI (двойная ось Y)
  p3 <- ggplot(all_results) +
    geom_line(aes(x = start_position, y = rf_distance, color = "RF расстояние"), linewidth = 1) +
    geom_line(aes(x = start_position, y = mci, 
                  color = "MCI (масштабировано)"), linewidth = 1) +
    scale_y_continuous(
      name = "RF расстояние",
      sec.axis = sec_axis(~ . * max(all_results$mci) / max(all_results$rf_distance),
                         name = "Mutual Clustering Info")
    ) +
    scale_color_manual(values = c("RF расстояние" = "steelblue", 
                                  "MCI (масштабировано)" = "darkgreen")) +
    labs(title = "Сравнение RF расстояния и Mutual Clustering Info",
         x = "Начальная позиция фрагмента",
         color = "Метрика") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  print(p3)
  ggsave("rf_mci_comparison.png", width = 10, height = 6, dpi = 300)
  
  # Сохраняем все результаты
  write.csv(all_results, "all_fragment_analysis_results_with_mci.csv", row.names = FALSE)
  
  # Выводим сводную статистику
  if (nrow(all_results) > 0) {
    cat("\nСводная статистика анализа:\n")
    
    cat("\nRF расстояние:\n")
    print(summary(all_results$rf_distance))
    
    cat("\nMutual Clustering Info:\n")
    print(summary(all_results$mci))
    
    # Находим оптимальные параметры
    # Для RF расстояния минимальное значение лучше
    optimal_rf <- all_results[which.min(all_results$rf_distance), ]
    
    # Для MCI максимальное значение лучше (больше сходства)
    optimal_mci <- all_results[which.max(all_results$mci), ]
    
    cat("\nОптимальные параметры по RF расстоянию (минимальное):\n")
    print(optimal_rf)
    
    cat("\nОптимальные параметры по MCI (максимальное):\n")
    print(optimal_mci)
    
    # Статистика по начальным позициям
    cat("\nОбщая статистика:\n")
    cat("Диапазон стартовых позиций:", min(all_results$start_position), "-", max(all_results$start_position), "\n")
    cat("Количество проанализированных позиций:", nrow(all_results), "\n")
  }
}