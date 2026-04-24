library(NMF)
library(Seurat)

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

intra_robust <- sapply(all_programs, function(prog) {
  same_sample <- all_programs[sapply(all_programs, function(x) x$sample == prog$sample)]
  overlaps <- sapply(same_sample, function(x) length(intersect(x$genes, prog$genes)))
  sum(overlaps >= 35) >= 2
})

inter_robust <- sapply(all_programs, function(prog) {
  other_samples <- all_programs[sapply(all_programs, function(x) x$sample != prog$sample)]
  any(sapply(other_samples, function(x) length(intersect(x$genes, prog$genes)) >= 10))
})

robust_programs <- all_programs[intra_robust & inter_robust]

program_similarity <- sapply(robust_programs, function(prog) {
  other_samples <- robust_programs[sapply(robust_programs, function(x) x$sample != prog$sample)]
  sum(sapply(other_samples, function(x) length(intersect(x$genes, prog$genes))))
})

sorted_programs <- robust_programs[order(program_similarity, decreasing = TRUE)]
unique_programs <- list()

for (prog_id in names(sorted_programs)) {
  prog <- sorted_programs[[prog_id]]
  
  is_redundant <- any(sapply(unique_programs, function(unique_prog) {
    length(intersect(unique_prog$genes, prog$genes)) > 10
  }))
  
  if (!is_redundant) {
    unique_programs[[prog_id]] <- prog
  }
}

cluster_mp <- function(programs, min_overlap = 10, max_genes = 200) {
  clusters <- list()
  remaining <- programs
  
  while (length(remaining) > 0) {
    connections <- sapply(remaining, function(p) {
      sum(sapply(remaining, function(q) {
        length(intersect(p$genes, q$genes)) >= min_overlap
      }))
    })
    
    seed_id <- names(which.max(connections))
    cluster <- list(seed_id)
    cluster_genes <- head(remaining[[seed_id]]$genes, 50) 

    added <- TRUE
    while (added && length(cluster_genes) < max_genes) {
      added <- FALSE
      candidates <- setdiff(names(remaining), unlist(cluster))
      
      best_match <- NULL
      best_ratio <- 0
      
      for (pid in candidates) {
        shared <- intersect(remaining[[pid]]$genes, cluster_genes)
        ratio <- length(shared)/length(remaining[[pid]]$genes)
        
        if (length(shared) >= min_overlap && ratio > best_ratio) {
          best_ratio <- ratio
          best_match <- pid
        }
      }
      
      if (!is.null(best_match)) {
        cluster <- c(cluster, best_match)
        new_genes <- setdiff(head(remaining[[best_match]]$genes, 30), cluster_genes)
        cluster_genes <- c(cluster_genes, new_genes)
        added <- TRUE
      }
    }
    
    clusters <- c(clusters, list(cluster))
    remaining <- remaining[setdiff(names(remaining), unlist(cluster))]
  }
  
  lapply(clusters, function(cluster_ids) {
    all_genes <- table(unlist(lapply(cluster_ids, function(id) {
      head(unique_programs[[id]]$genes, 30)
    })))
    
    sources <- do.call(rbind, lapply(cluster_ids, function(id) {
      p <- unique_programs[[id]]
      data.frame(
        program = id,
        sample = p$sample,
        cancer = p$cancer,
        stringsAsFactors = FALSE
      )
    }))
    
    list(
      programs = cluster_ids,
      sources = sources,
      genes = names(sort(all_genes, decreasing = TRUE)),
      gene_counts = as.vector(sort(all_genes, decreasing = TRUE)),
      top_genes = names(head(sort(all_genes, decreasing = TRUE), 50))
    )
  })
}

mps <- cluster_mp(unique_programs)
                          
filter_mps <- function(mps, multi_cancer = FALSE, min_prog_per_cancer = 1) {
  Filter(function(mp) {
    cancer_types <- unique(mp$sources$cancer)
    sample_contrib <- table(mp$sources$sample)
    cancer_contrib <- table(mp$sources$cancer)
    ribo_genes <- grep("^RPS|^RPL|^MRPS|^MRPL", mp$top_genes, value = TRUE)
    mito_genes <- grep("^MT-", mp$top_genes, value = TRUE)
    
    cond1 <- if (multi_cancer) length(cancer_types) > 1 else length(cancer_types) == 1
    cond2 <- length(mp$programs) >= 2
    cond3 <- if (multi_cancer) {
      min(cancer_contrib) >= min_prog_per_cancer
    } else {
      length(unique(mp$sources$sample[mp$sources$cancer == cancer_types])) >= 2
    }
    cond4 <- length(c(ribo_genes, mito_genes)) < 10
    
    passed <- cond1 && cond2 && cond3 && cond4
    passed
  }, mps)
}

filtered_mps_pan <- filter_mps(mps, multi_cancer = TRUE)
filtered_mps_single <- filter_mps(mps, multi_cancer = FALSE)

os_specific <- Filter(function(x) all(x$sources$cancer == "OS"), filtered_mps_single)
rb_specific <- Filter(function(x) all(x$sources$cancer == "RB"), filtered_mps_single)
pan_cancer <- Filter(function(x) {
  cancers <- unique(x$sources$cancer)
  "OS" %in% cancers && "RB" %in% cancers
}, filtered_mps_pan)

dir.create("data/malignant_programs", showWarnings = FALSE, recursive = TRUE)

saveRDS(list(
  all_mps = mps,
  os_specific = os_specific,
  rb_specific = rb_specific,
  pan_cancer = pan_cancer
), "data/malignant_programs/classified_MPs.rds")

generate_signature_matrix <- function(mp_list, prefix = "MP") {
  all_genes <- unique(unlist(lapply(mp_list, function(x) x$top_genes)))
  
  signature_matrix <- sapply(seq_along(mp_list), function(i) {
    mp <- mp_list[[i]]
    as.integer(all_genes %in% mp$top_genes)
  })
  
  colnames(signature_matrix) <- paste0(prefix, "_", seq_along(mp_list))
  rownames(signature_matrix) <- all_genes
  signature_matrix
}

if (length(os_specific) > 0) {
  write.csv(generate_signature_matrix(os_specific, "OS_MP"), "data/malignant_programs/OS_specific_MPs.csv")
}

if (length(rb_specific) > 0) {
  write.csv(generate_signature_matrix(rb_specific, "RB_MP"), "data/malignant_programs/RB_specific_MPs.csv")
}

if (length(pan_cancer) > 0) {
  write.csv(generate_signature_matrix(pan_cancer, "Pan_MP"), "data/malignant_programs/pan_cancer_MPs.csv")
}
