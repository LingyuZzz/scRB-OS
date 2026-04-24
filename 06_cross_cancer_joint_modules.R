library(Seurat)
library(dplyr)
library(ggplot2)
library(org.Hs.eg.db)

os_data <- readRDS('data/OS/os_MalignantOnly_816.rds')
rb_data <- readRDS('data/RB/rb_MalignantOnly_anno721.rds')

rb_data$cancer_type <- "RB"
os_data$cancer_type <- "OS"

combined <- merge(rb_data, y = os_data, add.cell.ids = c("RB", "OS"))

combined.list <- SplitObject(combined, split.by = "orig.ident")
combined.list <- lapply(combined.list, function(obj) {
  obj <- SCTransform(obj, verbose = FALSE)
  return(obj)
})

anchors <- FindIntegrationAnchors(object.list = combined.list)
combined <- IntegrateData(anchorset = anchors)

DefaultAssay(combined) <- "integrated"
combined <- ScaleData(combined)
combined <- RunPCA(combined, npcs = 30)
combined <- RunUMAP(combined, dims = 1:20)
combined <- FindNeighbors(combined, dims = 1:30)
combined <- FindClusters(combined, resolution = 0.5)

saveRDS(combined, 'data/OS+RB_mal_SCTcca_816.rds')

meta <- combined@meta.data
mix_cluster <- meta %>%
  group_by(seurat_clusters) %>%
  summarize(
    total = n(),
    rb = sum(cancer_type == "RB"),
    os = sum(cancer_type == "OS"),
    rb_ratio = rb / total,
    os_ratio = os / total
  ) %>%
  mutate(
    p_value = purrr::map2_dbl(rb, total, ~binom.test(.x, .y, p = 0.5)$p.value),
    dominant_type = ifelse(
      p_value >= 0.05, 
      "Mixed", 
      ifelse(rb_ratio > os_ratio, "RB_Dominant", "OS_Dominant")
    )
  ) %>% arrange(p_value)

write.csv(mix_cluster, 'results/mal_binom_test.csv', row.names = FALSE)

submodules_df <- read.csv("results/shared_submodules.csv")
sub1 <- submodules_df$gene[submodules_df$submodule == 1]
sub9 <- submodules_df$gene[submodules_df$submodule == 9]

combined <- AddModuleScore(combined, features = list(intersect(sub1, rownames(combined))), name = "ME_module")
combined <- AddModuleScore(combined, features = list(intersect(sub9, rownames(combined))), name = "RS_module")

rb_core_network <- c('RB1','CDK4','MDM2','CDK2','E2F1','E2F2','HDAC1','TFDP1','E2F4','CCND1','E2F3')  
rb_core_network <- mapIds(org.Hs.eg.db, keys = rb_core_network, column = "ENSEMBL", keytype = "SYMBOL", multiVals = "first")
rb_core_network <- unname(rb_core_network[!is.na(rb_core_network)])

exprData <- as.matrix(GetAssayData(combined, assay = "RNA", slot = "data"))

module1_score <- colMeans(exprData[intersect(sub1, rownames(exprData)), ])
module2_score <- colMeans(exprData[intersect(sub9, rownames(exprData)), ])
rb_score      <- colMeans(exprData[intersect(rb_core_network, rownames(exprData)), ])

cor_module1 <- cor.test(module1_score, rb_score, method = "spearman")
cor_module2 <- cor.test(module2_score, rb_score, method = "spearman")

cat("ME Module vs RB core spearman correlation:\n")
print(cor_module1)
cat("RS Module vs RB core spearman correlation:\n")
print(cor_module2)

df <- data.frame(
  RB_core = rb_score,
  Module1 = module1_score,
  Module2 = module2_score,
  Subtype = combined@meta.data$cancer_type
)

p1 <- ggplot(df, aes(x = RB_core, y = Module1, color = Subtype)) +
  geom_point(size = 0.3, alpha = 0.3) +  
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed", size = 0.8) +
  theme_minimal() + labs(x = "RB Core Activity", y = "ME Module Activity")

p2 <- ggplot(df, aes(x = RB_core, y = Module2, color = Subtype)) +
  geom_point(size = 0.3, alpha = 0.3) +  
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed", size = 0.8) +
  theme_minimal() + labs(x = "RB Core Activity", y = "RS Module Activity")

ggsave("plots/module_corr_scatter.pdf", plot = p1 + p2, width = 10, height = 5)
