if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
library(DESeq2)
library(readxl)
library(openxlsx)

counts_df <- read_excel("/Users/mortezaabyadeh/Desktop/gene_count.xlsx")
head(counts_df)

sample_names <- colnames(counts_df)[3:ncol(counts_df)] # Extract sample names from counts_df

# Create col_data with sample names and group information
col_data <- data.frame(
  sampleName = sample_names,
  group = factor(rep(c("CTRL", "Eto", "Dox", "FTY", "FTY_Dox", "FTY_Eto", "Eto_FTY", "Dox_FTY"), each = 3))
)

countData <- counts_df[, -c(1, 2)]

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

gene_names <- counts_df$gene_name


results_Eto_vs_CTRL <- cbind(gene_name = gene_names, results_Eto_vs_CTRL)
results_FTY_vs_CTRL <- cbind(gene_name = gene_names, results_FTY_vs_CTRL)
results_FTY_Eto_vs_CTRL <- cbind(gene_name = gene_names, results_FTY_Eto_vs_CTRL)
results_FTY_vs_Eto <- cbind(gene_name = gene_names, results_FTY_vs_Eto)
results_FTY_Eto_vs_Eto <- cbind(gene_name = gene_names, results_FTY_Eto_vs_Eto)
results_FTY_Eto_vs_FTY <- cbind(gene_name = gene_names, results_FTY_Eto_vs_FTY)

significant_eto_ctrl <- subset(results_Eto_vs_CTRL, padj < 0.05)
significant_fty_ctrl <- subset(results_FTY_vs_CTRL, padj < 0.05)
significant_fty_eto_ctrl <- subset(results_FTY_Eto_vs_CTRL, padj < 0.05)
significant_fty_eto <- subset(results_FTY_vs_Eto, padj < 0.05)
significant_fty_eto_eto <- subset(results_FTY_Eto_vs_Eto, padj < 0.05)
significant_fty_eto_fty <- subset(results_FTY_Eto_vs_FTY, padj < 0.05)





library(openxlsx)

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
