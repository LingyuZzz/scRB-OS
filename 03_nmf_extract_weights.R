library(NMF)
library(Seurat)
library(dplyr)
library(proxy)

final_nmf_results <- readRDS("data/NMF_rank_results/NMF_final_k3-9_results.rds")

sample_to_cancer <- c(
    'GSM5065165'='RB', 'GSM5065167'='RB', 'GSM7968797'='RB', 'GSM7968798'='RB',
    'GSM7968799'='RB', 'GSM7968800'='RB', 'GSM4952363'='OS', 'GSM4952364'='OS',
    'GSM4952365'='OS', 'GSM5155198'='OS', 'GSM5155199'='OS', 'GSM5155200'='OS',
    'GSM4600954'='OS', 'GSM4600966'='OS', 'GSM4600974'='OS', 'GSM4600978'='OS', 
    'GSM4600982'='OS', 'GSM4600986'='OS', 'GSM4600990'='OS', 'GSM4600998'='OS', 
    'GSM4601006'='OS', 'GSM4601014'='OS', 'GSM4601022'='OS', 'GSM7969979'='OS', 
    'GSM7969982'='OS'
)

extract_top_genes_with_weights <- function(nmf_res, n_genes = 50) {
  w_matrix <- nmf_res@fit@W
  lapply(seq_len(ncol(w_matrix)), function(i) {
    weights <- w_matrix[, i]
    names(weights) <- rownames(w_matrix)
    sorted_weights <- sort(weights, decreasing = TRUE)
    return(head(sorted_weights, n_genes))
  })
}

all_programs_data <- list()

for (sample_id in names(final_nmf_results)) {
  cancer_type <- sample_to_cancer[sample_id]
  if(is.na(cancer_type)) next
  
  sample_res <- final_nmf_results[[sample_id]]$nmf_results
  
  for (k_val in names(sample_res)) {
    top_genes_list <- extract_top_genes_with_weights(sample_res[[k_val]])
    
    for (prog_idx in seq_along(top_genes_list)) {
      gene_weights <- top_genes_list[[prog_idx]]
      prog_df <- data.frame(
        ProgramID = paste(sample_id, k_val, prog_idx, sep = "_"),
        Sample = sample_id,
        Cancer = cancer_type,
        K_Setting = k_val,
        Gene = names(gene_weights),
        Weight = as.numeric(gene_weights),
        Rank = 1:length(gene_weights),
        stringsAsFactors = FALSE
      )
      all_programs_data[[length(all_programs_data) + 1]] <- prog_df
    }
  }
}

all_programs <- all_programs_data
final_gene_table <- do.call(rbind, all_programs_data)

if(!dir.exists("data/malignant_programs")) dir.create("data/malignant_programs", recursive = TRUE)

write.csv(final_gene_table, "data/malignant_programs/NMF_All_Programs_Gene_Weights.csv", row.names = FALSE)

gene_lists <- lapply(all_programs, function(x) names(x$gene_weights))
prog_names <- names(all_programs)

calc_jaccard <- function(a, b) {
  length(intersect(a, b)) / length(union(a, b))
}

n <- length(gene_lists)
sim_mat <- matrix(0, nrow = n, ncol = n)
for (i in 1:(n-1)) {
  for (j in (i+1):n) {
    sim <- calc_jaccard(gene_lists[[i]], gene_lists[[j]])
    sim_mat[i, j] <- sim_mat[j, i] <- sim
  }
}
diag(sim_mat) <- 1
rownames(sim_mat) <- colnames(sim_mat) <- prog_names

dist_mat <- as.dist(1 - sim_mat)
hc <- hclust(dist_mat, method = "complete")
clusters <- cutree(hc, k = 10) 

mp_list <- list()

for (cl in unique(clusters)) {
  member_ids <- names(clusters[clusters == cl])
  member_info <- do.call(rbind, lapply(all_programs[member_ids], function(x) {
    data.frame(Sample = x$sample, Cancer = x$cancer, stringsAsFactors = FALSE)
  }))
  
  all_weights <- list()
  for (id in member_ids) {
    w <- all_programs[[id]]$gene_weights
    for (g in names(w)) {
      all_weights[[g]] <- c(all_weights[[g]], w[g])
    }
  }
  
  gene_summary <- data.frame(
    Gene = names(all_weights),
    Mean_Weight = sapply(all_weights, mean),
    Frequency = sapply(all_weights, length),
    stringsAsFactors = FALSE
  ) %>% arrange(desc(Mean_Weight))
  
  mp_list[[paste0("MP_", cl)]] <- list(
    genes = gene_summary,
    samples = unique(member_info$Sample),
    cancers = unique(member_info$Cancer),
    n_programs = length(member_ids)
  )
}

pan_cancer_mps <- Filter(function(x) {
  all(c("OS", "RB") %in% x$cancers) && x$n_programs >= 3
}, mp_list)

final_mp_table <- do.call(rbind, lapply(names(pan_cancer_mps), function(mp_name) {
  mp_data <- pan_cancer_mps[[mp_name]]
  df <- mp_data$genes
  df$MetaProgram = mp_name
  df$Contributing_Programs = mp_data$n_programs
  df$Cancers = paste(mp_data$cancers, collapse = ",")
  return(df)
}))

final_mp_table <- final_mp_table %>%
  group_by(MetaProgram) %>%
  slice_max(order_by = Mean_Weight, n = 100)

write.csv(final_mp_table, "data/malignant_programs/MetaProgram_Conserved_Weights.csv", row.names = FALSE)