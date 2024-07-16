if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
library(DESeq2)
library(readxl)
library(dplyr)
library(openxlsx)

counts_df <- read_excel("/Users/mortezaabyadeh/Desktop/gene_count.xlsx")
head(counts_df)
setwd("/Users/mortezaabyadeh/Desktop")
dim(counts_df)
sample_names <- colnames(counts_df)[3:ncol(counts_df)] # Extract sample names from counts_df
head(counts_df)
# Create col_data with sample names and group information
col_data <- data.frame(
  sampleName = sample_names,
  group = factor(rep(c("CTRL", "Eto", "Dox", "FTY", "FTY_Dox", "FTY_Eto", "Eto_FTY", "Dox_FTY"), each = 3))
)

countData <- counts_df[, -c(1, 2)]
rownames(countData) <- make.unique(counts_df[[1]])
rownames(countData)
head(countData)
dim(countData)
# Create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = col_data,
                              design = ~ group)
dds <- DESeq(dds)


results <- results(dds)


results_Eto_vs_CTRL <- results(dds, contrast = c("group", "Eto", "CTRL"))
results_FTY_vs_CTRL <- results(dds, contrast = c("group", "FTY", "CTRL"))
results_FTY_Eto_vs_CTRL <- results(dds, contrast = c("group", "FTY_Eto", "CTRL"))
results_FTY_vs_Eto <- results(dds, contrast = c("group", "FTY", "Eto"))
results_FTY_Eto_vs_Eto <- results(dds, contrast = c("group", "FTY_Eto", "Eto"))
results_FTY_Eto_vs_FTY <- results(dds, contrast = c("group", "FTY_Eto", "FTY"))

dim(counts_df)
head(counts_df)
ctrl_column <- counts_df[, c("HS_1", "HS_2", "HS_3")]
Eto_column <- counts_df[, c("HS_4", "HS_5", "HS_6")]
FTY_column <- counts_df[, c("HS_10", "HS_11", "HS_12")]
FTY_Eto_column <- counts_df[, c("HS_16", "HS_17", "HS_18")]

head(results_Eto_vs_CTRL)
head(ctrl_column)
head(Eto_column)
results_Eto_vs_CTRL_count <- cbind(results_Eto_vs_CTRL, ctrl_column, Eto_column)
head(results_Eto_vs_CTRL_count)

results_FTY_vs_CTRL_count <- cbind(results_FTY_vs_CTRL, ctrl_column, FTY_column)
results_FTY_Eto_vs_CTRL_count <- cbind(results_FTY_Eto_vs_CTRL, ctrl_column, FTY_Eto_column)
results_FTY_vs_Eto_count <- cbind(results_FTY_vs_Eto, Eto_column, FTY_column)
results_FTY_Eto_vs_Eto_count <- cbind(results_FTY_Eto_vs_Eto, Eto_column, FTY_Eto_column)
results_FTY_Eto_vs_FTY_count <- cbind(results_FTY_Eto_vs_FTY, FTY_column, FTY_Eto_column)

str(results_FTY_Eto_vs_FTY_count)

colnames(results_Eto_vs_CTRL_count)[7:12] <- c("CTRL-1", "CTRL-2", "CTRL-3", "Eto-1", "Eto-2", "Eto-3")
colnames(results_FTY_vs_CTRL_count)[7:12] <- c("CTRL-1", "CTRL-2", "CTRL-3", "FTY-1", "FTY-2", "FTY-3")
colnames(results_FTY_Eto_vs_CTRL_count)[7:12] <- c("CTRL-1", "CTRL-2", "CTRL-3", "FTY_Eto-1", "FTY_Eto-2", "FTY_Eto-3")
colnames(results_FTY_vs_Eto_count)[7:12] <- c("Eto-1", "Eto-2", "Eto-3","FTY-1","FTY-2","FTY-3")
colnames(results_FTY_Eto_vs_Eto_count)[7:12] <- c("Eto-1", "Eto-2", "Eto-3", "FTY_Eto-1", "FTY_Eto-2", "FTY_Eto-3")
colnames(results_FTY_Eto_vs_FTY_count)[7:12] <- c("FTY-1", "FTY-2", "FTY-3", "FTY_Eto-1", "FTY_Eto-2", "FTY_Eto-3")


results_Eto_vs_CTRL_count$Gene <- rownames(results_Eto_vs_CTRL_count)

results_Eto_vs_CTRL_df <- as.data.frame(results_Eto_vs_CTRL_count)

head(results_Eto_vs_CTRL_df)
write.xlsx(results_Eto_vs_CTRL_df, file = "results_Eto_vs_CTRL.xlsx")


results_FTY_vs_CTRL_count$Gene <- rownames(results_FTY_vs_CTRL_count)
results_FTY_Eto_vs_CTRL_count$Gene <- rownames(results_FTY_Eto_vs_CTRL_count)
results_FTY_vs_Eto_count$Gene <- rownames(results_FTY_vs_Eto_count)
results_FTY_Eto_vs_Eto_count$Gene <- rownames(results_FTY_Eto_vs_Eto_count)
results_FTY_Eto_vs_FTY_count$Gene <- rownames(results_FTY_Eto_vs_FTY_count)



results_FTY_vs_CTRL_df <- as.data.frame(results_FTY_vs_CTRL_count)
results_FTY_Eto_vs_CTRL_count_df <- as.data.frame(results_FTY_Eto_vs_CTRL_count)
results_FTY_vs_Eto_count_df <- as.data.frame(results_FTY_vs_Eto_count)
results_FTY_Eto_vs_Eto_count_df <- as.data.frame(results_FTY_Eto_vs_Eto_count)
results_FTY_Eto_vs_FTY_count_df <- as.data.frame(results_FTY_Eto_vs_FTY_count)

write.xlsx(results_FTY_vs_CTRL_df, file = "results_FTY_vs_CTRL.xlsx")
write.xlsx(results_FTY_Eto_vs_CTRL_count_df, file = "results_FTY_Eto_vs_CTRL.xlsx")
write.xlsx(results_FTY_vs_Eto_count_df, file = "results_FTY_vs_Eto.xlsx")
write.xlsx(results_FTY_Eto_vs_Eto_count_df, file = "results_FTY_Eto_vs_Eto.xlsx")
write.xlsx(results_FTY_Eto_vs_FTY_count_df, file = "results_FTY_Eto_vs_FTY.xlsx")


######################################################################################################### Let's do the next step in Data preprocessing script

groups <- rep(c("CTRL", "Eto", "Dox", "FTY", "FTY_Dox", "FTY_Eto", "Eto_FTY", "Dox_FTY"), each = 3)
replicate_numbers <- rep(1:3, times = length(groups) / 3)
unique_sample_names <- paste0(groups, "-", replicate_numbers)

colnames(counts_df)[3:ncol(counts_df)] <- unique_sample_names


head(counts_df)
relevant_counts_1 <- counts_df %>%
  select(gene_name = 1, starts_with("CTRL"), starts_with("Eto"))

relevant_counts_2 <- counts_df %>%
  select(gene_name = 1, starts_with("CTRL"), starts_with("FTY"))
relevant_counts_3 <- counts_df %>%
  select(gene_name = 1, starts_with("CTRL"), starts_with("FTY_Eto"))
relevant_counts_4 <- counts_df %>%
  select(gene_name = 1, starts_with("Eto"), starts_with("FTY"))
relevant_counts_5 <- counts_df %>%
  select(gene_name = 1, starts_with("FTY_Eto"), starts_with("Eto"))
relevant_counts_6 <- counts_df %>%
  select(gene_name = 1, starts_with("FTY_Eto"), starts_with("FTY"))





merged_results_Eto_vs_CTRL <- merge(results_Eto_vs_CTRL_df, relevant_counts_1, by = "gene_name")
merged_results_FTY_vs_CTRL <- merge(results_Eto_vs_CTRL_df, relevant_counts_2, by = "gene_name")
merged_results_FTY_Eto_vs_CTRL <- merge(results_Eto_vs_CTRL_df, relevant_counts_3, by = "gene_name")
merged_results_FTY_Eto <- merge(results_Eto_vs_CTRL_df, relevant_counts_4, by = "gene_name")
merged_results_FTY_Eto_vs_Eto <- merge(results_Eto_vs_CTRL_df, relevant_counts_5, by = "gene_name")
merged_results_FTY_Eto_vs_FTY <- merge(results_Eto_vs_CTRL_df, relevant_counts_6, by = "gene_name")


write.xlsx(merged_results_Eto_vs_CTRL, file = "merged_results_Eto_vs_CTRL.xlsx")




significant_eto_ctrl <- subset(results_Eto_vs_CTRL, padj < 0.05)
significant_fty_ctrl <- subset(results_FTY_vs_CTRL, padj < 0.05)
significant_fty_eto_ctrl <- subset(results_FTY_Eto_vs_CTRL, padj < 0.05)
significant_fty_eto <- subset(results_FTY_vs_Eto, padj < 0.05)
significant_fty_eto_eto <- subset(results_FTY_Eto_vs_Eto, padj < 0.05)
significant_fty_eto_fty <- subset(results_FTY_Eto_vs_FTY, padj < 0.05)






excel_file <- "/Users/mortezaabyadeh/Desktop/DEGs.xlsx"

wb <- createWorkbook()

addWorksheet(wb, sheetName = "Eto_vs_CTRL")
writeData(wb, sheet = "Eto_vs_CTRL", significant_eto_ctrl)

addWorksheet(wb, sheetName = "FTY_vs_CTRL")
writeData(wb, sheet = "FTY_vs_CTRL", significant_fty_ctrl)

addWorksheet(wb, sheetName = "FTY_Eto_vs_CTRL")
writeData(wb, sheet = "FTY_Eto_vs_CTRL", significant_fty_eto_ctrl)

addWorksheet(wb, sheetName = "FTY_vs_Eto")
writeData(wb, sheet = "FTY_vs_Eto", significant_fty_eto)

addWorksheet(wb, sheetName = "FTY_Eto_vs_Eto")
writeData(wb, sheet = "FTY_Eto_vs_Eto", significant_fty_eto_eto)

addWorksheet(wb, sheetName = "FTY_Eto_vs_FTY")
writeData(wb, sheet = "FTY_Eto_vs_FTY", significant_fty_eto_fty)

saveWorkbook(wb, excel_file)



filterAndWrite <- function(data, sheet_name) {
  Down <- subset(data, log2FoldChange < 0)
  
  Up <- subset(data, log2FoldChange > 0)
  
  addWorksheet(wb, sheetName = paste0(sheet_name, "_Down"))
  writeData(wb, sheet = paste0(sheet_name, "_Down"), Down)
  
  addWorksheet(wb, sheetName = paste0(sheet_name, "_Up"))
  writeData(wb, sheet = paste0(sheet_name, "_Up"), Up)
}

filterAndWrite(significant_eto_ctrl, "Eto_vs_CTRL")
filterAndWrite(significant_fty_ctrl, "FTY_vs_CTRL")
filterAndWrite(significant_fty_eto_ctrl, "FTY_Eto_vs_CTRL")
filterAndWrite(significant_fty_eto, "FTY_vs_Eto")
filterAndWrite(significant_fty_eto_eto, "FTY_Eto_vs_Eto")
filterAndWrite(significant_fty_eto_fty, "FTY_Eto_vs_FTY")

saveWorkbook(wb, "/Users/mortezaabyadeh/Desktop/DEGs_final.xlsx")





### proteome and transcriptome comparision




