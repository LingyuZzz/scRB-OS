SAMPLE_TO_CANCER <- c(
  'GSM5065165' = 'RB', 'GSM5065167' = 'RB', 'GSM7968797' = 'RB',
  'GSM7968798' = 'RB', 'GSM7968799' = 'RB', 'GSM7968800' = 'RB',
  'GSM4952363' = 'OS', 'GSM4952364' = 'OS', 'GSM4952365' = 'OS',
  'GSM5155198' = 'OS', 'GSM5155199' = 'OS', 'GSM5155200' = 'OS',
  'GSM4600954' = 'OS', 'GSM4600966' = 'OS', 'GSM4600974' = 'OS',
  'GSM4600978' = 'OS', 'GSM4600982' = 'OS', 'GSM4600986' = 'OS',
  'GSM4600990' = 'OS', 'GSM4600998' = 'OS', 'GSM4601006' = 'OS',
  'GSM4601014' = 'OS', 'GSM4601022' = 'OS', 'GSM7969979' = 'OS',
  'GSM7969982' = 'OS'
)

N_RB_SAMPLES <- sum(SAMPLE_TO_CANCER == "RB")
N_OS_SAMPLES <- sum(SAMPLE_TO_CANCER == "OS")

RB_MARKER_GENES <- list(
  Rods = c("RHO", "PDE6A", "CNGA1", "GNAT1", "GNB1", "SAG", "GNGT1", "NR2E3", "PDE6G"),
  Cones = c("ARR3", "GNGT2", "GUCA1C", "GNAT2"),
  Microglia = c("HLA-DPA1", "HLA-DPB1", "HLA-DRA", "C1QA", "AIF1"),
  Astrocyte = c("GFAP", "AQP4", "CD44"),
  Muller_glia = c("GLUL", "APOE", "SPP1"),
  Fibroblast = c("COL3A1", "LUM", "DCN", "FBLN1"),
  TILs = c("IL7R", "CD3D", "NKG7"),
  Retinoma_like = c("CDCA7"),
  Proliferating_PhrD = c("TOP2A", "UBE2C", "BIRC5", "TPX2")
)

OS_MARKER_GENES <- list(
  Osteoclast = c("CTSK", "MMP9"),
  MSCs = c("CXCL12", "SFRP2", "MME"),
  Pericyte = c("ACTA2", "RGS5"),
  Fibroblast = c("COL1A1", "LUM", "COL3A1"),
  Myeloid = c("CD74", "CD14", "FCGR3A", "LYZ", "S100A9", "C1QA", "APOE"),
  Osteoblast = c("ALPL", "RUNX2", "CLEC11A"),
  TILs = c("IL7R", "CD3D", "NKG7"),
  Endothelial = c("VWF", "PECAM1", "CAV1", "CLDN5", "EGFL7")
)

DATA_PATHS <- list(
  os_raw = "data/OS/OS_addPhase_424.rds",
  rb_raw = "data/RB/RB_annotated_final.rds",
  rb_mal = "data/RB/rb_MalignantOnly_anno721.rds",
  os_mal = "data/OS/os_MalignantOnly_clustered_817.rds",
  nmf_results = "data/NMF_rank_results/NMF_final_k3-9_results.rds",
  infercnv_rb_results = "RB_inferCNV_results/",
  os_cellchat = "data/OS/os_cellchat_final_nboot.rds",
  rb_cellchat = "data/RB/rb_cellchat_final_nboot.rds",
  gene_pos = "data/gencode_v47_gene_pos.infercnv.txt"
)

CORES_TO_USE <- 32
CELLCHAT_NBOOT <- 100
SCENIC_CORES <- 32

RB_REF_CELLS <- c("Rods", "Cones", "TILs", "Astrocyte", "Microglia", "Muller glia")
OS_REF_CELLS <- c("MSCs", "Pericyte", "Fibroblast", "TILs", "Endothelial")

NMF_K_RANGE <- 3:9
NMF_NRUN <- 10
NMF_METHOD <- "brunet"

ROBUST_OVERLAP_THRESHOLD <- 35
INTER_ROBUST_THRESHOLD <- 10
REDUNDANCY_THRESHOLD <- 10
SIMILARITY_THRESHOLD_LOW <- 20
SIMILARITY_THRESHOLD_HIGH <- 70

SURVIAL_MIN_SAMPLES <- 20
SURVIVAL_EXPRESSION_CUTOFF <- 0.5
