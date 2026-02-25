setwd("/Users/mortezaabyadeh/Documents/npc-fty/transcriptome npc cells/DNA break transcriptome with different filteration")


library(readxl)


library(data.table)

fpkm <- fread("gene_fpkm.xls", header = TRUE, data.table = FALSE)

dim(fpkm)

gene <- read_excel("gene.xlsx")
head(gene)
head(fpkm)

fpkm1 <- fpkm[, c("gene_id", "S1", "S2", "S3", "S7", "S8", "S9", "gene_name", "gene_biotype")]
dim(fpkm1)
head(gene)
fpkm2 <- fpkm1[fpkm1$gene_id %in% gene$ensembl_id, ]
fpkm3 <- fpkm1[fpkm1$gene_name %in% gene$gene_name, ]

head(fpkm2)
dim(gene)
dim(fpkm2)
head(fpkm3)
dim(fpkm3)

#library(openxlsx)

### only protein coding
fpkm2_coding <- fpkm2[fpkm2$gene_biotype == "protein_coding", ]
dim(fpkm2_coding)
write.xlsx(fpkm2_coding, "fpkm2coding.xlsx", rowNames = TRUE)


write.xlsx(fpkm2, "fpkm2.xlsx", rowNames = TRUE)
write.xlsx(fpkm3, "fpkm3.xlsx", rowNames = TRUE)



######### DEG analysis #############
setwd("/Users/mortezaabyadeh/Desktop")
counts <- fread("gene_count.xls", header = TRUE, data.table = FALSE)

# Set gene_id as rownames
rownames(counts) <- counts$gene_id

# -----------------------------
# Keep only expression columns + gene_name
# -----------------------------
counts_sub <- counts[, c("gene_id",
                         paste0("S",1:12),
                         "gene_name")]

# Set rownames again (important after subsetting)
rownames(counts_sub) <- counts_sub$gene_id

# Create gene_name vector
gene_name_vec <- counts_sub$gene_name
names(gene_name_vec) <- counts_sub$gene_id

# Keep only count matrix
counts_only <- counts_sub %>%
  dplyr::select(paste0("S",1:12))

rownames(counts_only) <- counts_sub$gene_id

# -----------------------------
# Metadata
# -----------------------------
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

# -----------------------------
# DESeq2
# -----------------------------
dds <- DESeqDataSetFromMatrix(
  countData = counts_only,
  colData = metadata,
  design = ~ group
)

dds$group <- relevel(dds$group, ref = "CTRL")

# Light filtering
keep <- rowSums(counts(dds) >= 10) >= 2
dds <- dds[keep, ]

cat("Number of genes after filtering:", nrow(dds), "\n")

dds <- DESeq(dds)

# -----------------------------
# Pairwise comparisons
# -----------------------------
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
  
  # normalized counts
  norm_counts <- counts(dds, normalized = TRUE) %>%
    as.data.frame() %>%
    rownames_to_column("ensembl_id")
  
  df <- res_df %>%
    left_join(norm_counts, by = "ensembl_id")
  
  # samples per group
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

saveWorkbook(wb, "DE_all_pairwise_clean-new1.xlsx", overwrite = TRUE)


###### FPKM only for DEGs

gene <- read_excel("fpkm2.xlsx")
deg <- read_excel("DEGs.xlsx")
dim(gene)
dim(deg)
degfpkm <- gene[gene$gene_id %in% deg$ensembl_id, ]
dim(degfpkm)
write.xlsx(degfpkm, "degfpkm.xlsx", rowNames = TRUE)
