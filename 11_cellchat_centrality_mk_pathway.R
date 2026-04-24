library(CellChat)
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)

os_cellchat <- readRDS('data/OS/os_cellchat_final_nboot.rds')
rb_cellchat <- readRDS('data/RB/rb_cellchat_final_nboot.rds')

os_cellchat <- netAnalysis_computeCentrality(os_cellchat, slot.name = "netP")
rb_cellchat <- netAnalysis_computeCentrality(rb_cellchat, slot.name = "netP")

focus_signals <- c('ANGPTL', 'APP', 'CADM', 'CD46', 'CDH', 'JAM', 'LAMININ', 'MIF', 'MK', 'MPZ', 'NCAM', 'NOTCH', 'PTN', 'SEMA4', 'THBS', 'VEGF')

get_pathway_centrality <- function(cellchat_obj, pathways, dataset_name) {
  df <- data.frame(pathway = pathways, outdeg = NA, betweenness = NA)
  for (p in pathways) {
    if (p %in% names(cellchat_obj@netP$centr)) {
      p_data <- cellchat_obj@netP$centr[[p]]
      df$outdeg[df$pathway == p] <- sum(p_data$outdeg, na.rm = TRUE)
      df$betweenness[df$pathway == p] <- sum(p_data$betweenness, na.rm = TRUE)
    }
  }
  df$degree_rank <- rank(-df$outdeg, na.last = "keep")
  df$betweenness_rank <- rank(-df$betweenness, na.last = "keep")
  df$avg_rank <- (df$degree_rank + df$betweenness_rank) / 2
  colnames(df)[-1] <- paste0(colnames(df)[-1], "_", dataset_name)
  return(df)
}

os_centr <- get_pathway_centrality(os_cellchat, focus_signals, "OS")
rb_centr <- get_pathway_centrality(rb_cellchat, focus_signals, "RB")

combined_ranks <- merge(os_centr, rb_centr, by = "pathway")
combined_ranks$mean_avg_rank <- (combined_ranks$avg_rank_OS + combined_ranks$avg_rank_RB) / 2
conservative_hubs <- combined_ranks[order(combined_ranks$mean_avg_rank), ]

write.csv(conservative_hubs, 'results/Joint_conservative_hubs.csv', row.names = FALSE)

p1 <- ggplot(conservative_hubs, aes(x = degree_rank_OS, y = degree_rank_RB)) +
  geom_point(aes(color = mean_avg_rank, size = 1/mean_avg_rank), alpha = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  geom_label_repel(data = conservative_hubs[1:7, ], aes(label = pathway, fill = mean_avg_rank),
                   color = "white", fontface = "bold") +
  scale_color_viridis_c(option = "plasma", direction = -1) +
  scale_fill_viridis_c(option = "plasma", direction = -1, guide = "none") +
  scale_size_continuous(guide = "none") + theme_minimal() + coord_fixed(ratio = 1) +
  labs(title = "Comparison of Pathway Degree Centrality", x = "OS Degree Rank", y = "RB Degree Rank")

p2 <- ggplot(conservative_hubs, aes(x = betweenness_rank_OS, y = betweenness_rank_RB)) +
  geom_point(aes(color = mean_avg_rank, size = 1/mean_avg_rank), alpha = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  geom_label_repel(data = conservative_hubs[1:7, ], aes(label = pathway, fill = mean_avg_rank),
                   color = "white", fontface = "bold") +
  scale_color_viridis_c(option = "plasma", direction = -1) +
  scale_fill_viridis_c(option = "plasma", direction = -1, guide = "none") +
  scale_size_continuous(guide = "none") + theme_minimal() + coord_fixed(ratio = 1) +
  labs(title = "Comparison of Pathway Betweenness Centrality", x = "OS Betweenness Rank", y = "RB Betweenness Rank")

library(patchwork)
ggsave('plots/Joint_signal_rank.pdf', p1 + p2, width = 12, height = 6)

os_group.cellType <- c("mal", "non_mal", "non_mal", "mal", "non_mal", "non_mal", 'non_mal', 'non_mal', 'mal', 'non_mal')
names(os_group.cellType) <- levels(os_cellchat@idents)
rb_group.cellType <- c("non_mal", "non_mal", "mal", "non_mal", "non_mal", "mal", 'non_mal', 'mal', 'non_mal', 'non_mal', 'mal', 'non_mal')
names(rb_group.cellType) <- levels(rb_cellchat@idents)

get_common_lr_data <- function(pathway_name, os_obj, rb_obj) {
  os_lr_df <- subsetCommunication(os_obj, signaling = pathway_name)
  rb_lr_df <- subsetCommunication(rb_obj, signaling = pathway_name)
  os_lr_df <- os_lr_df[os_group.cellType[os_lr_df$source] == "mal" | os_group.cellType[os_lr_df$target] == "mal", ]
  rb_lr_df <- rb_lr_df[rb_group.cellType[rb_lr_df$source] == "mal" | rb_group.cellType[rb_lr_df$target] == "mal", ]
  
  common_interactions <- intersect(os_lr_df$interaction_name, rb_lr_df$interaction_name)
  if (length(common_interactions) == 0) return(NULL)
  
  result_df <- data.frame(interaction_name = common_interactions, os_prob = numeric(length(common_interactions)), rb_prob = numeric(length(common_interactions)))
  for (i in seq_along(common_interactions)) {
    lr_name <- common_interactions[i]
    result_df$os_prob[i] <- max(os_lr_df$prob[os_lr_df$interaction_name == lr_name], na.rm = TRUE)
    result_df$rb_prob[i] <- max(rb_lr_df$prob[rb_lr_df$interaction_name == lr_name], na.rm = TRUE)
  }
  return(result_df)
}

all_lr_data <- lapply(focus_signals, function(p) {
  res <- get_common_lr_data(p, os_cellchat, rb_cellchat)
  if (!is.null(res)) res$pathway <- p
  res
})
combined_data <- do.call(rbind, all_lr_data)

if (!is.null(combined_data)) {
  plot_data <- combined_data %>%
    pivot_longer(cols = c(os_prob, rb_prob), names_to = "dataset", values_to = "probability") %>%
    mutate(dataset = ifelse(dataset == "os_prob", "OS", "RB"),
           lr_pathway = paste0(interaction_name, " (", pathway, ")"))

  p_bar <- ggplot(plot_data, aes(x = lr_pathway, y = 1/log10(probability + 1e-5), fill = dataset)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.75) +
    scale_fill_manual(values = c("OS" = "#ff9d9a", "RB" = "#1ABC9C"), name = "Cancer Type") +
    labs(x = "Ligand-Receptor Pair", y = "Communication Probability (1/log10)", title = "Shared Ligand-Receptor Pairs") +
    theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + coord_flip()
    
  ggsave("plots/Joint_LR_Pair.pdf", p_bar, height = 10, width = 7, dpi = 300)
}

targets_mal_os <- names(os_group.cellType)[os_group.cellType == "mal"]
sources_non_mal_os <- setdiff(names(os_group.cellType), targets_mal_os)

pdf("plots/os_MK_chord.pdf", width = 8, height = 8) 
netVisual_individual(os_cellchat, signaling = "MK", pairLR.use = "MDK_NCL", layout = "chord", 
                     group = os_group.cellType, sources.use = sources_non_mal_os, targets.use = targets_mal_os)
dev.off()

targets_mal_rb <- names(rb_group.cellType)[rb_group.cellType == "mal"]
sources_non_mal_rb <- setdiff(names(rb_group.cellType), targets_mal_rb)

pdf("plots/rb_MK_chord.pdf", width = 8, height = 8) 
netVisual_individual(rb_cellchat, signaling = "MK", pairLR.use = "MDK_NCL", layout = "chord", 
                     group = rb_group.cellType, sources.use = sources_non_mal_rb, targets.use = targets_mal_rb)
dev.off()
