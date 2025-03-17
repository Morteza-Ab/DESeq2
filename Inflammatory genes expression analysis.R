library(readxl)
library(dplyr)
setwd("/Users/mortezaabyadeh/Documents/npc-fty/senescence omics data/Data vis/data to analyze")
file1 <- "clean_transcriptome_Eto_vs_CTRL.xlsx"
file2 <- "clean_transcriptome_FTY_Eto_vs_Eto.xlsx"

df1 <- read_excel(file1, sheet = "DEGs") %>% select(1:7)
df2 <- read_excel(file2, sheet = "DEGs") %>% select(1:7)
head(df1)
head(df2)

dim(df1)
dim(df2)
#colnames(df1) <- c("Gene", paste0("Eto_vs_CTRL_", colnames(df1)[-1]))
#colnames(df2) <- c("Gene", paste0("FTY_Eto_vs_Eto_", colnames(df2)[-1]))


common_genes <- intersect(df1$Gene, df2$Gene)

df1_common <- df1 %>% filter(Gene %in% common_genes)
df2_common <- df2 %>% filter(Gene %in% common_genes)

dim(df1_common)
dim(df2_common)

merged_df <- merge(df1_common, df2_common, by = "Gene")
head(merged_df)

merged_df <- merged_df %>% select(-c(`Eto-1.y`, `Eto-2.y`, `Eto-3.y`))
head(merged_df)


biomarker <- read_excel("/Users/mortezaabyadeh/Documents/npc-fty/senescence omics data/Data vis/inflammatory genes cyto/Inflammatory biomarker measurment.xlsx")
head(biomarker)

common_genes1 <- intersect(merged_df$Gene, biomarker$Gene)

filtered_merge <- merged_df %>% filter(Gene %in% common_genes1)

dim(merged_df)
dim(filtered_merge)

filtered_merge <- filtered_merge %>%
  rename(`FTY_Eto-1` = `FTY-Eto-1`, 
         `FTY_Eto-2` = `FTY-Eto-2`, 
         `FTY_Eto-3` = `FTY-Eto-3`,
         `Eto-1` = `Eto-1.x`, 
         `Eto-2` = `Eto-2.x`, 
         `Eto-3` = `Eto-3.x`)


head(filtered_merge)

######################### Add the fingolimod #########################
df3 <- df1 <- read_excel("clean_transcriptome_FTY_vs_CTRL.xlsx", sheet = "DEGs") %>% select(1:7)

head(df3)

common <- intersect(df3$Gene, filtered_merge$Gene)
df3_common <- df3 %>% filter(Gene %in% common)
merged_common <- filtered_merge %>% filter(Gene %in% common)

merged_df1 <- merge(merged_common, df3_common, by = "Gene")

head(merged_df1)

merged_df1 <- merged_df1 %>% select(-c(`CTRL-1.y`, `CTRL-2.y`, `CTRL-3.y`))

dim(merged_df1)


merged_df1 <- merged_df1 %>%
  rename(`FTY_Eto-1` = `FTY-Eto-1`, 
         `FTY_Eto-2` = `FTY-Eto-2`, 
         `FTY_Eto-3` = `FTY-Eto-3`)

merged_df <- merged_df %>%
  rename(`FTY_Eto-1` = `FTY-Eto-1`, 
         `FTY_Eto-2` = `FTY-Eto-2`, 
         `FTY_Eto-3` = `FTY-Eto-3`)



merged_df <- merged_df %>%
  rename(`Eto-1` = `Eto_1`, 
         `Eto-2` = `Eto_2`, 
         `Eto-3` = `Eto_3`)

colnames(merged_df)

library(pheatmap)


heatmap_maker <- function(df) {
  df_filtered <- df[complete.cases(df), ]
  
  # Set the Gene column as row names and remove it from the dataset
  rownames(df_filtered) <- df_filtered$Gene
  df_filtered <- df_filtered[, -1]
  
  i_matrix <- as.matrix(df_filtered)
  
  i_log <- log2(i_matrix + 1)
  
  i_scaled <- t(scale(t(i_log), scale = TRUE))
  
  group_names <- sub("-.*", "", colnames(i_scaled))
  
  i_annotation_col <- data.frame(Group = group_names)
  rownames(i_annotation_col) <- colnames(i_scaled)
  
  pheatmap(i_scaled,
           cluster_cols = TRUE,
           cluster_rows = TRUE,
           scale = "row",
           color = colorRampPalette(c("green", "black", "red"))(100),
           main = "Heatmap of Normalized Counts",
           fontsize = 13,
           fontsize_row = 12,
           fontsize_col = 13,
           border_color = NA,
           annotation_col = i_annotation_col,
           show_rownames = TRUE,
           width = 10,
           height = 7
  )
}


heatmap_maker(merged_df1)
heatmap_maker(filtered_merge)
dim(merged_df)
