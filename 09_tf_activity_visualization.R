library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)

os <- readRDS('data/OS/os_MalignantOnly_clustered_817.rds')
rb <- readRDS('data/RB/rb_MalignantOnly_anno721.rds')

os_activity <- read.csv('data/scenic/os_tf_activity_matrix.csv', row.names = 1)
rb_activity <- read.csv('data/scenic/rb_tf_activity_matrix.csv', row.names = 1)

colnames(os_activity) <- trimws(sub("\\.\\.\\..*$", "", colnames(os_activity)))
colnames(rb_activity) <- trimws(sub("\\.\\.\\..*$", "", colnames(rb_activity)))

common_cells_os <- intersect(rownames(os@meta.data), rownames(os_activity))
os <- subset(os, cells = common_cells_os)

common_cells_rb <- intersect(rownames(rb@meta.data), rownames(rb_activity))
rb <- subset(rb, cells = common_cells_rb)

make_tf_scores_long <- function(activity, meta, tf_group1, tf_group2, tf_single) {
  act_df <- as.data.frame(activity, stringsAsFactors = FALSE)
  act_df$cell <- rownames(act_df)
  
  match_cols_ci <- function(requested, cols) {
    cols_up <- toupper(cols)
    sapply(requested, function(r) {
      i <- which(cols_up == toupper(r))
      if (length(i) >= 1) cols[i[1]] else NA_character_
    }, USE.NAMES = FALSE)
  }
  
  present1 <- na.omit(match_cols_ci(tf_group1, colnames(act_df)))
  present2 <- na.omit(match_cols_ci(tf_group2, colnames(act_df)))
  present_single <- match_cols_ci(tf_single, colnames(act_df))
  
  act_df$Group1 <- if(length(present1) > 0) rowMeans(act_df[, present1, drop = FALSE], na.rm = TRUE) else NA_real_
  act_df$Group2 <- if(length(present2) > 0) rowMeans(act_df[, present2, drop = FALSE], na.rm = TRUE) else NA_real_
  
  for (i in seq_along(tf_single)) {
    tf <- tf_single[i]; actual <- present_single[i]
    if (!is.na(actual) && actual != tf) act_df[[tf]] <- act_df[[actual]]
  }
  
  sel_cols <- c("cell", tf_single, "Group1", "Group2")
  scores_df <- act_df[, intersect(sel_cols, colnames(act_df)), drop = FALSE]
  
  meta_df <- as.data.frame(meta, stringsAsFactors = FALSE)
  meta_df$cell <- rownames(meta_df)
  
  merged <- merge(scores_df, meta_df, by = "cell", all.x = TRUE)
  
  long <- pivot_longer(merged, cols = any_of(c(tf_single, "Group1", "Group2")), 
                       names_to = "TF", values_to = "Activity") %>%
    mutate(TF_group = case_when(
      TF == "Group1" ~ "Group1",
      TF == "Group2" ~ "Group2",
      TF %in% tf_single ~ "Core",
      TRUE ~ "Other"
    ))
  return(long)
}

tf_group1 <- c('EZH2','BRCA1','E2F8')
tf_group2 <- c('TFDP1','E2F2','E2F1')
tf_single <- c('RAD21','E2F7')

os_long <- make_tf_scores_long(os_activity, os@meta.data, tf_group1, tf_group2, tf_single)
rb_long <- make_tf_scores_long(rb_activity, rb@meta.data, tf_group1, tf_group2, tf_single)

p_os <- ggplot(os_long, aes(x = TF, y = Activity, fill = TF_group)) +
  geom_violin(trim = FALSE, scale = "width") +
  geom_boxplot(width = 0.1, outlier.size = 0.3, alpha = 0.5) +
  facet_wrap(~subtype, scales = "free_y") +
  theme_bw(base_size = 12) + coord_flip() +
  labs(title = "OS TF/Group Activity by Subtype")

p_rb <- ggplot(rb_long, aes(x = TF, y = Activity, fill = TF_group)) +
  geom_violin(trim = FALSE, scale = "width") +
  geom_boxplot(width = 0.1, outlier.size = 0.3, alpha = 0.5) +
  facet_wrap(~subtype, scales = "free_y") +
  theme_bw(base_size = 12) + coord_flip() +
  labs(title = "RB TF/Group Activity by Subtype")

ggsave("plots/os_mal_group_activity_vln.pdf", p_os, width = 10, height = 5, units = "in")
ggsave("plots/rb_mal_group_activity_vln.pdf", p_rb, width = 10, height = 5, units = "in")
