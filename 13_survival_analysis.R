library(survival)
library(survminer)
library(dplyr)
library(org.Hs.eg.db)

expr_data <- read.delim("data/TCGA/TARGET-OS.star_fpkm-uq.tsv.gz", check.names = FALSE)
surv_data <- read.delim("data/TCGA/TARGET-OS.survival.tsv.gz")

expr_data$Ensembl_ID <- sub("\\..*", "", expr_data$Ensembl_ID)

surv_samples <- surv_data$sample
expr_samples <- colnames(expr_data)[-1] 
common_samples <- intersect(expr_samples, surv_samples)
expr_data <- expr_data[, c("Ensembl_ID", common_samples)]

analyze_single_gene_survival <- function(gene_symbol) {
  gene_ensembl <- unname(mapIds(org.Hs.eg.db, keys = gene_symbol, column = "ENSEMBL", keytype = "SYMBOL", multiVals = "first"))
  if (is.na(gene_ensembl)) return(NULL)

  gene_row <- expr_data[expr_data$Ensembl_ID == gene_ensembl, ]
  if (nrow(gene_row) == 0) return(NULL)

  expr_df <- data.frame(
    sample_id = colnames(gene_row)[-1],
    expression = as.numeric(gene_row[1, -1])
  )

  merged_data <- surv_data %>%
    dplyr::rename(sample_id = sample) %>%
    inner_join(expr_df, by = "sample_id") %>%
    filter(!is.na(expression) & !is.na(OS.time) & OS.time > 0)

  if (nrow(merged_data) < 20) return(NULL)

  median_expr <- median(merged_data$expression, na.rm = TRUE)
  merged_data$expression_group <- ifelse(merged_data$expression > median_expr, "High", "Low")
  
  surv_obj <- Surv(time = merged_data$OS.time, event = merged_data$OS)
  fit <- survfit(surv_obj ~ expression_group, data = merged_data)

  km_plot <- ggsurvplot(
    fit, data = merged_data,
    pval = TRUE, conf.int = TRUE, risk.table = TRUE,
    legend.labs = c("High Expression", "Low Expression"),
    legend.title = gene_symbol,
    title = paste("Survival Analysis for", gene_symbol),
    xlab = "Time (days)", ylab = "Overall Survival Probability",
    palette = c("#E7B800", "#2E9FDF")
  )
  
  ggsave(paste0("plots/", gene_symbol, "_survival_curve.pdf"), km_plot$plot, width = 6, height = 6)
  
  return(km_plot)
}

analyze_single_gene_survival("NCL")
