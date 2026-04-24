library(monocle3)
library(Seurat)
library(org.Hs.eg.db)
library(dplyr)
library(tidyr)  
library(ggplot2)

rb_data <- readRDS('data/RB/rb_MalignantOnly_anno721.rds')
os_data <- readRDS('data/OS/os_MalignantOnly_clustered_817.rds')

rb_data <- RenameIdents(
    object = rb_data,
    'subtype1' = 'Mitotic subtype',
    'subtype2' = 'Stress-adaptive subtype',
    'subtype3' = 'DNA damage response subtype',
    'subtype4' = 'Replication stress subtype'
)
rb_data@meta.data$subtype <- Idents(rb_data)

expr_mat <- GetAssayData(rb_data, layer = "counts")  
cell_meta <- rb_data@meta.data
gene_meta <- data.frame(gene_short_name = rownames(expr_mat), row.names = rownames(expr_mat))
rb_cds <- new_cell_data_set(expr_mat, cell_metadata = cell_meta, gene_metadata = gene_meta)
rb_cds <- preprocess_cds(rb_cds, num_dim = 50)   
reducedDims(rb_cds)$PCA <- rb_data@reductions$pca@cell.embeddings
rb_cds <- cluster_cells(rb_cds, reduction_method = "PCA")
reducedDims(rb_cds)$UMAP <- rb_data@reductions$umap@cell.embeddings
rb_cds <- cluster_cells(rb_cds, reduction_method = "UMAP")
rb_cds@clusters$UMAP$clusters <- rb_data$subtype
rb_cds <- learn_graph(rb_cds)
rb_cds <- order_cells(rb_cds, root_cells = colnames(rb_cds)[rb_cds$subtype == "Mitotic subtype"])

p_rb <- plot_cells(rb_cds, color_cells_by = "pseudotime", label_groups_by_cluster = TRUE,
                   label_leaves = TRUE, label_branch_points = TRUE) +
                   theme(aspect.ratio = 1) + coord_fixed(ratio = 1)
ggsave(filename = 'plots/pdf/rb_pseudotime.pdf', plot = p_rb, width = 6, height = 6, dpi = 300, bg = 'white')

expr_mat <- GetAssayData(os_data, layer = "counts")  
cell_meta <- os_data@meta.data
gene_meta <- data.frame(gene_short_name = rownames(expr_mat), row.names = rownames(expr_mat))
os_cds <- new_cell_data_set(expr_mat, cell_metadata = cell_meta, gene_metadata = gene_meta)
os_cds <- preprocess_cds(os_cds, num_dim = 50)   
reducedDims(os_cds)$PCA <- os_data@reductions$pca@cell.embeddings
os_cds <- cluster_cells(os_cds, reduction_method = "PCA")
reducedDims(os_cds)$UMAP <- os_data@reductions$umap@cell.embeddings
os_cds <- cluster_cells(os_cds, reduction_method = "UMAP")
os_cds@clusters$UMAP$clusters <- os_data$subtype
os_cds <- learn_graph(os_cds)
os_cds <- order_cells(os_cds, root_cells = colnames(os_cds)[os_cds$subtype == "Mitotic subtype"])

p_os <- plot_cells(os_cds, raster = T, color_cells_by = "pseudotime", label_groups_by_cluster = TRUE,
                   label_leaves = TRUE, label_branch_points = TRUE) +
                   theme(aspect.ratio = 1) + coord_fixed(ratio = 1)
ggsave(filename = 'plots/pdf/os_pseudotime.pdf', plot = p_os, width = 6, height = 6, dpi = 300, bg = 'white')

pdata <- as.data.frame(colData(rb_cds))
pdata$cell_id <- rownames(pdata)
pdata$pseudotime <- pseudotime(rb_cds) 
colnames(pdata)[colnames(pdata) == "Proliferation1"] <- "Proliferation"
colnames(pdata)[colnames(pdata) == "DNA_Repair1"] <- "DNA_Repair"
colnames(pdata)[colnames(pdata) == "Glycolysis1"] <- "Glycolysis"
colnames(pdata)[colnames(pdata) == "Apoptosis1"] <- "Apoptosis"
colnames(pdata)[colnames(pdata) == "P53_PATHWAY1"] <- "P53_PATHWAY"
colnames(pdata)[colnames(pdata) == "HYPOXIA1"] <- "HYPOXIA"
score_cols <- c("Proliferation", "DNA_Repair", "Glycolysis", "Apoptosis", 'P53_PATHWAY', 'HYPOXIA')
df_long <- pdata %>% select(cell_id, pseudotime, all_of(score_cols)) %>% pivot_longer(cols = all_of(score_cols), names_to = "hallmark", values_to = "score")

prb <- ggplot(df_long, aes(x = pseudotime, y = score)) +
  geom_smooth(method = "loess", span = 0.5, se = FALSE, linewidth = 1.3, color = "#1B7837") +
  facet_wrap(~hallmark, scales = "free_y", ncol = 3) +
  labs(x = "Pseudotime", y = "Module score") +
  theme_minimal(base_size = 13) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        strip.background = element_rect(fill = "#F7F7F7", color = NA),
        strip.text = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 13), panel.spacing = unit(1.3, "lines"))
ggsave(filename = 'plots/rb_hallmark_pseudotime.pdf', plot = prb, width = 10, height = 4, dpi = 300, bg = 'white')

submodules <- read.csv("results/shared_submodules.csv")
sub1 <- submodules$gene[submodules$submodule == 1]
sub9 <- submodules$gene[submodules$submodule == 9]

DefaultAssay(rb_data) <- "RNA"
rb_data <- AddModuleScore(object = rb_data, features = list(intersect(sub9, rownames(rb_data))), name = "RS_module", ctrl = 80)
rb_data <- AddModuleScore(object = rb_data, features = list(intersect(sub1, rownames(rb_data))), name = "ME_module", ctrl = 80)

DefaultAssay(os_data) <- "RNA"
os_data <- AddModuleScore(object = os_data, features = list(intersect(sub9, rownames(os_data))), name = "RS_module", ctrl = 80)
os_data <- AddModuleScore(object = os_data, features = list(intersect(sub1, rownames(os_data))), name = "ME_module", ctrl = 80)

pdata_rb <- as.data.frame(colData(rb_cds))
pdata_rb$cell_id <- rownames(pdata_rb)
pdata_rb$pseudotime <- pseudotime(rb_cds) 
pdata_rb$ME_module <- rb_data@meta.data$ME_module1
pdata_rb$RS_module <- rb_data@meta.data$RS_module1
df_long_rb <- pdata_rb %>% select(cell_id, pseudotime, ME_module, RS_module) %>% pivot_longer(cols = c('ME_module', 'RS_module'), names_to = "Module", values_to = "score")

p1 <- ggplot(df_long_rb, aes(x = pseudotime, y = score, color = Module)) +
  geom_smooth(method = "loess", se = F) + theme_minimal() +
  labs(x = "RB Pseudotime", y = "Expression Score", color = "Module") + theme(legend.position = "right")

pdata_os <- as.data.frame(colData(os_cds))
pdata_os$cell_id <- rownames(pdata_os)
pdata_os$pseudotime <- pseudotime(os_cds) 
pdata_os$ME_module <- os_data@meta.data$ME_module1
pdata_os$RS_module <- os_data@meta.data$RS_module1
df_long_os <- pdata_os %>% select(cell_id, pseudotime, ME_module, RS_module) %>% pivot_longer(cols = c('ME_module', 'RS_module'), names_to = "Module", values_to = "score")

p2 <- ggplot(df_long_os, aes(x = pseudotime, y = score, color = Module)) +
  geom_smooth(method = "loess", se = F) + theme_minimal() +
  labs(x = "OS Pseudotime", y = "Expression Score", color = "Module") + theme(legend.position = "right")

ps_all <- (p1 | p2) 
ggsave("plots/modules_score_pseudotime.pdf", plot = ps_all, width = 10, height = 4, dpi = 300)
