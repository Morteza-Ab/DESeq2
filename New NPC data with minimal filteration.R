setwd("/Users/mortezaabyadeh/Desktop")

library(DESeq2)
library(dplyr)
library(openxlsx)
library(EnhancedVolcano)


counts <- read.table(
  "gene_count.xls",
  header = TRUE,
  row.names = 1,
  sep = "\t",
  check.names = FALSE
)

counts_sub <- counts[, c(1:12, which(colnames(counts) == "gene_name"))]

gene_name_vec <- counts_sub$gene_name
names(gene_name_vec) <- rownames(counts_sub)

counts_only <- counts_sub %>% dplyr::select(S1:S12)
rownames(counts_only) <- rownames(counts_sub)



metadata <- data.frame(
  sample = colnames(counts_only),
  group = factor(c(
    rep("CTRL", 3),
    rep("CTRL_T", 3),
    rep("NPC", 3),
    rep("NPC_T", 3)
  ))
)

rownames(metadata) <- metadata$sample
metadata$sample <- NULL


dds <- DESeqDataSetFromMatrix(
  countData = counts_only,
  colData = metadata,
  design = ~ group
)

dds$group <- relevel(dds$group, ref = "CTRL")


keep <- rowSums(counts(dds) >= 1) >= 2
dds <- dds[keep, ]

cat("Number of genes after filtering:", nrow(dds), "\n")


dds <- DESeq(dds)


groups <- levels(dds$group)
pairwise <- combn(groups, 2, simplify = FALSE)

wb <- createWorkbook()

for(x in pairwise){
  
  comp_name <- paste(x[1], "vs", x[2])
  
  res <- results(
    dds,
    contrast = c("group", x[1], x[2]),
    independentFiltering = FALSE
  )
  
  res_df <- as.data.frame(res) %>%
    rownames_to_column("ensembl_id") %>%
    mutate(gene_name = gene_name_vec[ensembl_id])
  
# normalized count
  norm_counts <- counts(dds, normalized = TRUE) %>%
    as.data.frame() %>%
    rownames_to_column("ensembl_id")
  
  df <- res_df %>%
    left_join(norm_counts, by = "ensembl_id")
  

  samples_grp1 <- rownames(metadata)[metadata$group == x[1]]
  samples_grp2 <- rownames(metadata)[metadata$group == x[2]]
  

  geom_mean <- function(x) exp(mean(log(x + 1)))
  
  df$geomean_grp1 <- apply(df[, samples_grp1], 1, geom_mean)
  df$geomean_grp2 <- apply(df[, samples_grp2], 1, geom_mean)
  
  df$fold_change <- df$geomean_grp2 / df$geomean_grp1
  df$log2FoldChange_calc <- log2(df$fold_change)
  
  df_out <- df %>%
    dplyr::select(
      ensembl_id,
      gene_name,
      all_of(colnames(norm_counts)[-1]),
      geomean_grp1,
      geomean_grp2,
      fold_change,
      log2FoldChange_calc,
      log2FoldChange,
      pvalue,
      padj
    )
  
  addWorksheet(wb, sheetName = comp_name)
  writeData(wb, sheet = comp_name, df_out)
}

saveWorkbook(wb, "DE_all_pairwise_clean.xlsx", overwrite = TRUE)
