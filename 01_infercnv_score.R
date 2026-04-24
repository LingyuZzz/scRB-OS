.libPaths("/mnt/public5/pancancercol/miniconda3/envs/R441/lib/R/library")
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(reticulate)
})

# ============================================================================
# Step 1: Compute CNV scores using infercnvpy (Python)
# ============================================================================
py_run_string("
import infercnvpy as cnv
import scipy.sparse as sp
import scanpy as sc

def process_infercnv(h5ad_path, output_csv):
    print(f'Reading {h5ad_path} ...')
    adata = sc.read_h5ad(h5ad_path)
    
    if 'X_cnv' in adata.obsm:
        X_cnv = adata.obsm['X_cnv']
        if sp.issparse(X_cnv):
            adata.obsm['X_cnv'] = X_cnv.toarray()
            
    if 'cnv_score' not in adata.obs.columns:
        print('Computing cnv_score...')
        cnv.tl.cnv_score(adata)
        
    out = adata.obs[['cnv_score']].copy()
    out.index.name = 'barcode'
    out.to_csv(output_csv)

process_infercnv('OS_inferCNV_results/infercnv_working_result.h5ad', 'OS_inferCNV_results/cnv_score.csv')
process_infercnv('RB_inferCNV_results/infercnv_working_result.h5ad', 'RB_inferCNV_results/cnv_score.csv')
")

# ============================================================================
# Step 2: Add CNV scores to Seurat objects and visualize
# ============================================================================
os_all <- readRDS('data/OS/OS_AnnotateMalignant_624.rds')
rb_all <- readRDS('data/RB/RB_AnnotateMalignant_624.rds')

os_cnv_df <- read.csv("OS_inferCNV_results/cnv_score.csv", row.names = 1)
rb_cnv_df <- read.csv("RB_inferCNV_results/cnv_score.csv", row.names = 1)

os_all$cnv_score <- os_cnv_df[colnames(os_all), "cnv_score"]
rb_all$cnv_score <- rb_cnv_df[colnames(rb_all), "cnv_score"]

os_all$cnv_score_plot <- ifelse(os_all$celltypes == "Osteoclast", NA, os_all$cnv_score)

os_plot <- FeaturePlot(
  os_all, features = "cnv_score", reduction = "umap", 
  min.cutoff = "q10", max.cutoff = "q90", order = TRUE, raster = TRUE
) + scale_colour_gradientn(colours = c("grey85", "#fd8d3c", "#bd0026"), na.value = "grey90") +
  theme(aspect.ratio = 1) + coord_fixed(ratio = 1)

rb_plot <- FeaturePlot(
  rb_all, features = "cnv_score", reduction = "umap", 
  min.cutoff = "q10", max.cutoff = "q95", order = TRUE, raster = TRUE
) + scale_colour_gradientn(colours = c("grey85", "#fd8d3c", "#bd0026"), na.value = "grey90") +
  theme(aspect.ratio = 1) + coord_fixed(ratio = 1)

ggsave('plots/os_cnv_score.pdf', os_plot)
ggsave('plots/rb_cnv_score.pdf', rb_plot)

