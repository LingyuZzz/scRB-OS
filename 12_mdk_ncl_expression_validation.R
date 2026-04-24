library(Seurat)
library(dplyr)
library(ggplot2)
library(org.Hs.eg.db)
library(patchwork)
library(viridis)
library(ggpubr)
library(RColorBrewer)
library(tidyr)
library(pheatmap)
library(AnnotationDbi)

os_mal <- readRDS('data/OS/os_MalignantOnly_clustered_817.rds')
rb_mal <- readRDS('data/RB/rb_MalignantOnly_anno721.rds')

os_mal$ME_module <- os_mal$ME_module1
os_mal$RS_module <- os_mal$RS_module1
rb_mal$ME_module <- rb_mal$ME_module1
rb_mal$RS_module <- rb_mal$RS_module1

ncl_ensembl <- "ENSG00000115053"

plot_ncl_correlation <- function(seu, name_prefix) {
  seu$NCL_expr <- FetchData(seu, vars = ncl_ensembl)[,1]
  df <- data.frame(RS_module = seu$RS_module, NCL_expr = seu$NCL_expr)
  df <- df[df$NCL_expr > 0, ]
  
  cor_res <- cor.test(df$RS_module, df$NCL_expr, method = "spearman")
  rho <- round(cor_res$estimate, 3)
  pval <- signif(cor_res$p.value, 3)
  
  set.seed(123)
  df_sample <- if (nrow(df) > 2000) df[sample(1:nrow(df), 2000), ] else df
  
  p <- ggplot(df_sample, aes(x = RS_module, y = NCL_expr)) +
    geom_point(shape = 21, color = "#3B4CC0", fill = "#A1D99B", size = 1) +
    geom_smooth(method = "lm", color = "#E64B35", fill = "#FDBF6F", se = TRUE, linetype = "solid") +
    xlab("DR module score") + ylab("NCL expression") +
    theme_bw() + theme(plot.title = element_text(hjust = 0.5), panel.grid.minor = element_blank()) +
    annotate("text", x = min(df_sample$RS_module), y = max(df_sample$NCL_expr), 
             label = paste0("rho = ", rho, "\n", "p = ", pval), hjust = 0, size = 4, color = "black")
             
  ggsave(paste0("plots/DR_vs_NCL_", name_prefix, ".pdf"), p, width = 5, height = 5, dpi = 300)
}

plot_ncl_correlation(os_mal, "OS")
plot_ncl_correlation(rb_mal, "RB")

rb_reg <- read.csv("data/rb_ncl_regulators_filtered.csv")
os_reg <- read.csv("data/os_ncl_regulators_filtered.csv")

os_tfs <- unname(mapIds(org.Hs.eg.db, keys = os_reg$TF, keytype = "SYMBOL", column = "ENSEMBL", multiVals = "first"))
rb_tfs <- unname(mapIds(org.Hs.eg.db, keys = rb_reg$TF, keytype = "SYMBOL", column = "ENSEMBL", multiVals = "first"))

os_tfs <- os_tfs[!is.na(os_tfs) & os_tfs %in% rownames(os_mal)]
rb_tfs <- rb_tfs[!is.na(rb_tfs) & rb_tfs %in% rownames(rb_mal)]

plot_tf_violin <- function(seu, tfs, name_prefix) {
  ensembl_to_symbol <- mapIds(org.Hs.eg.db, keys = tfs, keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first")
  subtypes <- unique(seu$subtype)
  plot_list <- list()
  
  for (i in seq_along(subtypes)) {
    target_subtype <- subtypes[i]
    cells_in_subtype <- WhichCells(seu, expression = subtype == target_subtype)
    sub_obj <- subset(seu, cells = cells_in_subtype)
    expr_mat <- as.data.frame(t(as.matrix(sub_obj@assays$RNA@data[tfs, ])))
    expr_mat$cell <- rownames(expr_mat)
    
    expr_long <- expr_mat %>% pivot_longer(cols = all_of(tfs), names_to = "TF", values_to = "Expression") %>%
      mutate(SYMBOL = ensembl_to_symbol[TF])
      
    p <- ggplot(expr_long, aes(x = SYMBOL, y = Expression, fill = SYMBOL)) +
      geom_violin(trim = FALSE, scale = "width") + theme_classic() +
      ggtitle(target_subtype) + xlab("") + ylab("Expression") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none", plot.title = element_text(hjust = 0.5))
    plot_list[[i]] <- p
  }
  
  combined_plot <- wrap_plots(plot_list, ncol = 1, guides = "collect")
  ggsave(paste0("plots/", name_prefix, "_ncl_TF.pdf"), combined_plot, width = 10, height = 2 * length(subtypes), dpi = 300, bg = "white")
}

plot_tf_violin(os_mal, os_tfs, "os")
plot_tf_violin(rb_mal, rb_tfs, "rb")

plot_tf_heatmap <- function(seu, tfs, name_prefix) {
  ensembl_to_symbol <- mapIds(org.Hs.eg.db, keys = tfs, keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first")
  avg_expr <- AverageExpression(seu, features = tfs, group.by = "subtype", slot = "data")$RNA
  rownames(avg_expr) <- ensembl_to_symbol[rownames(avg_expr)]
  
  pheatmap(t(avg_expr), cluster_rows = FALSE, cluster_cols = TRUE, border_color = NA,
           fontsize_row = 10, fontsize_col = 8, scale = "column", cellwidth = 15, cellheight = 15,
           treeheight_col = 10, angle_col = 45, legend_position = "left",
           filename = paste0("plots/", name_prefix, "_ncl_tf_heatmap.pdf"), width = 10, height = 6)
}

plot_tf_heatmap(os_mal, os_tfs, "os")
plot_tf_heatmap(rb_mal, rb_tfs, "rb")
