library(NMF)
library(Seurat)
library(pheatmap)
library(viridis)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)

final_nmf_results <- readRDS("data/NMF_rank_results/NMF_final_k3-9_results.rds")
all_samples <- names(final_nmf_results)

sample_to_cancer <- c(
    'GSM5065165'='RB', 'GSM5065167'='RB', 'GSM7968797'='RB', 'GSM7968798'='RB',
    'GSM7968799'='RB', 'GSM7968800'='RB', 'GSM4952363'='OS', 'GSM4952364'='OS', 
    'GSM4952365'='OS', 'GSM5155198'='OS', 'GSM5155199'='OS', 'GSM5155200'='OS', 
    'GSM4600954'='OS', 'GSM4600966'='OS', 'GSM4600974'='OS', 'GSM4600978'='OS', 
    'GSM4600982'='OS', 'GSM4600986'='OS', 'GSM4600990'='OS', 'GSM4600998'='OS', 
    'GSM4601006'='OS', 'GSM4601014'='OS', 'GSM4601022'='OS', 'GSM7969979'='OS', 
    'GSM7969982'='OS'
)

missing_samples <- setdiff(all_samples, names(sample_to_cancer))
if (length(missing_samples) > 0) {
  stop(paste("Missing mappings:", paste(missing_samples, collapse = ", ")))
}

extract_top_genes <- function(nmf_res, n_genes = 50) {
  lapply(seq(ncol(nmf_res@fit@W)), function(i) {
    coef <- nmf_res@fit@W[, i]
    names(coef) <- rownames(nmf_res@fit@W)
    head(names(sort(coef, decreasing = TRUE)), n_genes)
  })
}

all_programs <- list()
for (sample_id in names(final_nmf_results)) {
  cancer_type <- sample_to_cancer[sample_id]
  sample_res <- final_nmf_results[[sample_id]]$nmf_results
  
  for (k in names(sample_res)) {
    top_genes <- extract_top_genes(sample_res[[k]])
    
    for (i in seq_along(top_genes)) {
      program_id <- paste(sample_id, k, i, sep = "_")
      all_programs[[program_id]] <- list(
        sample = sample_id,
        cancer = cancer_type,
        k = k,
        genes = top_genes[[i]]
      )
    }
  }
}

prog_names <- names(all_programs)
gene_list <- lapply(all_programs, function(x) x$genes)
meta_info <- data.frame(
  ProgramID = prog_names,
  Sample = sapply(all_programs, function(x) x$sample),
  Cancer = sapply(all_programs, function(x) x$cancer),
  stringsAsFactors = FALSE
)

calc_jaccard <- function(list1, list2) {
  inter <- length(intersect(list1, list2))
  uni <- length(union(list1, list2))
  return(inter / uni)
}

n <- length(gene_list)
jaccard_mat <- matrix(0, nrow = n, ncol = n)
rownames(jaccard_mat) <- prog_names
colnames(jaccard_mat) <- prog_names

for (i in 1:n) {
  for (j in i:n) {
    val <- calc_jaccard(gene_list[[i]], gene_list[[j]])
    jaccard_mat[i, j] <- val
    jaccard_mat[j, i] <- val
  }
}

cancer_colors <- c("OS" = "#E41A1C", "RB" = "#377EB8")
samples <- unique(meta_info$Sample)
sample_colors <- setNames(scales::hue_pal()(length(samples)), samples)

ha <- rowAnnotation(
  Cancer = meta_info$Cancer,
  Sample = meta_info$Sample,
  col = list(Cancer = cancer_colors, Sample = sample_colors),
  show_annotation_name = TRUE
)

pdf("data/malignant_programs/all_Programs_Jaccard_Heatmap.pdf", width = 10, height = 8)

magma_cols <- rev(magma(100))
col_fun = colorRamp2(
  c(0, seq(0.01, 0.5, length.out = 100)), 
  c("white", magma_cols)
)

Heatmap(jaccard_mat, 
        name = "Jaccard Index",
        column_title = "Similarity between NMF Programs",
        col = col_fun,
        left_annotation = ha,
        cluster_rows = TRUE, 
        cluster_columns = TRUE,
        clustering_distance_rows = function(m) as.dist(1 - m),
        clustering_distance_columns = function(m) as.dist(1 - m),
        clustering_method_rows = "complete",
        show_row_names = FALSE, 
        show_column_names = FALSE,
        border = TRUE)

dev.off()
