library(NMF)
library(Seurat)

cores_to_use <- 32
options(mc.cores = cores_to_use)
nmf.options("cores" = cores_to_use)

combined <- readRDS('data/OS+RB_mal_SCTcca_816.rds')

output_dir <- "data/NMF_rank_results"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

result_file <- file.path(output_dir, "NMF_final_k3-9_results.rds")
temp_result_file <- file.path(output_dir, "NMF_temp_results.rds")

if (file.exists(temp_result_file)) {
  message("Loading existing temporary results...")
  final_nmf_results <- readRDS(temp_result_file)
} else if (file.exists(result_file)) {
  message("Loading existing final results...")
  final_nmf_results <- readRDS(result_file)
} else {
  final_nmf_results <- list()
}

all_samples <- unique(combined$orig.ident)
processed_samples <- names(final_nmf_results)
remaining_samples <- setdiff(all_samples, processed_samples)

message("Already processed: ", length(processed_samples), " samples")
message("Remaining samples to process: ", length(remaining_samples))

for (sample_id in remaining_samples) {
  message("\n", paste(rep("=", 50), collapse = ""))
  message("Processing sample: ", sample_id, " (", 
          which(remaining_samples == sample_id), "/", 
          length(remaining_samples), ")")
  message(paste(rep("=", 50), collapse = ""))
  
  start_time <- Sys.time()
  
  tryCatch({
    seu_sub <- subset(combined, subset = orig.ident == sample_id)
    
    if(ncol(seu_sub) < 50) {
      message("  Skipping sample ", sample_id, ": too few cells (", ncol(seu_sub), ")")
      final_nmf_results[[sample_id]] <- list(
        status = "skipped",
        reason = paste("too few cells:", ncol(seu_sub))
      )
      saveRDS(final_nmf_results, temp_result_file)
      next
    }
    
    DefaultAssay(seu_sub) <- "RNA"
    seu_sub <- NormalizeData(seu_sub)
    seu_sub <- FindVariableFeatures(seu_sub, nfeatures = 2000)
    
    var_features <- VariableFeatures(seu_sub)
    if(length(var_features) < 10) {
      message("  Skipping sample ", sample_id, ": insufficient variable features (", length(var_features), ")")
      final_nmf_results[[sample_id]] <- list(
        status = "skipped",
        reason = paste("insufficient variable features:", length(var_features))
      )
      saveRDS(final_nmf_results, temp_result_file)
      next
    }
    
    expr_mat <- as.matrix(GetAssayData(seu_sub, slot = "data")[var_features, ])
    expr_mat <- t(scale(t(expr_mat)))
    expr_mat[expr_mat < 0] <- 0
    
    nmf_res <- list()
    for(k in 3:9) {
      message("  Running NMF with k=", k)
      nmf_res[[as.character(k)]] <- nmf(expr_mat, rank = k, method = "brunet", nrun = 10)
    }
    
    final_nmf_results[[sample_id]] <- list(
      orig.ident = sample_id,
      nmf_results = nmf_res,
      cell_barcodes = colnames(seu_sub),
      status = "completed",
      processing_time = as.numeric(difftime(Sys.time(), start_time, units = "mins")),
      n_cells = ncol(seu_sub),
      n_features = length(var_features)
    )
    
    message("  Successfully processed sample: ", sample_id, 
            " (Time: ", round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 1), " mins)")
    
    saveRDS(final_nmf_results, temp_result_file)
    message("  Temporary results saved for sample: ", sample_id)
    
  }, error = function(e) {
    message("  ERROR processing sample ", sample_id, ": ", e$message)
    final_nmf_results[[sample_id]] <- list(
      status = "error",
      error_message = e$message,
      sample_id = sample_id
    )
    saveRDS(final_nmf_results, temp_result_file)
    message("  Error results saved for sample: ", sample_id)
  }, warning = function(w) {
    message("  WARNING in sample ", sample_id, ": ", w$message)
  })
}

if (file.exists(temp_result_file)) {
  file.rename(temp_result_file, result_file)
  message("Final results saved to: ", result_file)
} else {
  saveRDS(final_nmf_results, result_file)
  message("Final results saved to: ", result_file)
}

message("Processing completed. ", 
      sum(sapply(final_nmf_results, function(x) x$status == "completed")), 
      " samples processed successfully.")

final_nmf_results <- readRDS("data/NMF_rank_results/NMF_final_k3-9_results.rds")
all_samples <- names(final_nmf_results)