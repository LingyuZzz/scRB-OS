library(CellChat)
library(Seurat)
library(future)
options(future.globals.maxSize = 10 * 1024^3) # 10GB
library(org.Hs.eg.db)

os_data <- readRDS('data/OS/OS_addPhase_424.rds')
rb_data <- readRDS('data/RB/RB_annotated_final.rds')
rb_mal <- readRDS('data/RB/rb_MalignantOnly_anno721.rds')
os_mal <- readRDS('data/OS/os_MalignantOnly_clustered_817.rds')

rb_mal@meta.data$subtype <- Idents(rb_mal)
rb_data$subtype <- as.character(rb_data$celltypes)
rb_data$subtype[rb_data$celltypes == "Proliferating PhrD"] <- NA
rb_data$subtype[colnames(rb_mal)] <- as.character(rb_mal@meta.data$subtype)

os_data$subtype <- as.character(os_data$celltypes)
os_data$subtype[os_data$celltypes == "Osteoblast"] <- NA
os_data$subtype[colnames(os_mal)] <- as.character(os_mal@meta.data$subtype)

run_cellchat <- function(seurat_obj, out_path) {
  cells_use <- rownames(seurat_obj@meta.data)[!is.na(seurat_obj$subtype)]
  sub_obj <- subset(seurat_obj, cells = cells_use)
  
  symbol_ids <- mapIds(org.Hs.eg.db, keys = rownames(sub_obj), keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first")
  keep_genes <- !is.na(symbol_ids)
  sub_obj <- sub_obj[keep_genes, ]
  
  data.input <- sub_obj[["RNA"]]@data
  rownames(data.input) <- unname(symbol_ids[keep_genes])
  
  meta <- data.frame(labels = sub_obj$subtype, row.names = colnames(sub_obj)) 
  cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
  cellchat@DB <- CellChatDB.human
  cellchat <- subsetData(cellchat)
  
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat, type = "triMean", nboot = 100)
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  
  saveRDS(cellchat, out_path)
  return(cellchat)
}

os_cellchat <- run_cellchat(os_data, "data/OS/os_cellchat_final_nboot.rds")
rb_cellchat <- run_cellchat(rb_data, "data/RB/rb_cellchat_final_nboot.rds")
