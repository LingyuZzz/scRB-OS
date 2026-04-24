Source code for the cross-cancer single-cell analysis of Osteosarcoma and Retinoblastoma.

This repository contains code of analyses in the manuscript. The scripts are arranged sequentially from 00 to 14.

**Pre_processing_and_CNV**
`00_setup_config.R` configures the environment and loads reference parameters. `01_infercnv_score.R` calculates inferred copy number variations (inferCNV).

**Meta_Programs_Identification**
These scripts (`02_run_nmf.R` to `05_mp_jaccard_similarity.R`) generate and identify robust intra-tumoral meta-programs from NMF programs. NMF is executed on malignant subsets, and Jaccard similarity is used to map these programs across OS and RB, defining cross-cancer conserved programs.

**Joint_Analysis_and_Trajectories**
`06_cross_cancer_joint_modules.R` identifies the shared Mitotic Execution and Replication Stress Adaptation modules across both cancers. `07_developmental_trajectories.R` infers pseudotime differentiation trajectories underlying these identified functional states using Monocle3.

**Regulatory_and_Communication_Networks**
`08_scenic_tf_enrichment.py` and `09_tf_activity_visualization.R` extract and visualize core transcription factor regulons via SCENIC. The CellChat scripts (`10_build_cellchat_network.R` to `12_mdk_ncl_expression_validation.R`) construct interaction networks and compute communication centrality, identifying the highly conserved MK (Midkine) signaling pathway and validating the MDK-NCL ligand-receptor interaction.

**Survival_Analysis**
`13_survival_analysis.R` conducts Kaplan-Meier and Cox proportional hazards survival analyses for key targets using TCGA data.
