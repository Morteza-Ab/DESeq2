install.packages("WGCNA", dependencies = TRUE)
setwd("/Users/mortezaabyadeh/Desktop")
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")

BiocManager::install("preprocessCore")
library(WGCNA)
BiocManager::install("GEOquery")

install.packages(c(
  "BiocManager",      # For Bioconductor package management
  "ggplot2",          # For advanced visualization
  "reshape2",         # For data manipulation
  "igraph",           # For network analysis and visualization
  "pheatmap",         # For heatmap visualization
  "RColorBrewer",     # For color palettes
  "corrplot",         # For correlation plots
  "ggrepel"           # For non-overlapping labels in plots
))
BiocManager::install(c(
  "limma",            # For normalization and QC
  "DESeq2",           # For data transformation
  "clusterProfiler",  # For functional enrichment analysis
  "org.Hs.eg.db",     # Human gene annotations
  "GO.db",            # Gene Ontology database
  "STRINGdb"          # For protein-protein interaction networks
))
library(WGCNA)
library(ggplot2)
library(ggrepel)
library(reshape2)
library(pheatmap)
library(igraph)
library(clusterProfiler)
library(org.Hs.eg.db)
library(STRINGdb)
library(DESeq2)

enableWGCNAThreads(nThreads = 8)
BiocManager::install("GEOquery")
library(GEOquery) 
### data normalization ###

getGEOSuppFiles("GSE261875")

# List all files in the downloaded folder
files <- list.files("GSE261875", full.names = TRUE)

# Read the count data file (first file in the list)
count_data <- read.delim(files[1])

# Extract sample information (metadata)
sample_info <- getGEO("GSE261875", GSEMatrix = TRUE)
sample_info <- pData(sample_info[[1]])


sample_info <- data.frame(
  SampleID = sample_info$description,
  Treatment = gsub("\\d+$", "", sample_info$description)  # Remove trailing numbers
)
rownames(sample_info) <- sample_info[[1]]

# Save sample_info table
write.table(sample_info, "sample_info.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# Match sample IDs between sample_info and count_data
# Keep gene ID columns and matching samples
count_data <- cbind(count_data[, c(1, 2)], count_data[, rownames(sample_info)])

# Verify that sample IDs are matched correctly
identical(rownames(sample_info), colnames(count_data)[-c(1, 2)])

#-----------------------------------------------
# STEP 3: Prepare count matrix
#-----------------------------------------------

# Remove rows with missing values
count_data <- na.omit(count_data)

# Remove duplicate gene entries (keep first occurrence)
count_data <- count_data[!duplicated(count_data[[2]]), ]

# Set gene IDs as row names
rownames(count_data) <- count_data[[2]]

# Create a separate gene annotation dataframe
gene_data <- count_data[, c(2, 1)]
rownames(gene_data) <- gene_data[[1]]

# Extract only the count columns (remove gene ID columns)
count_data <- count_data[, -c(1, 2)]

# Filter low-count genes
# Keep genes with >5 counts in at least 50% of samples
keep <- rowSums(count_data > 5) >= 0.5 * ncol(count_data)
count_data <- count_data[keep, ]

# Update gene annotations to match filtered counts
gene_data <- gene_data[rownames(count_data), ]

# Save gene data
write.table(gene_data, "gene_data.tsv", sep = "\t", row.names = FALSE, quote = FALSE)



### STEP 4: Normalize count data using DESeq2
#-----------------------------------------------

# Create DESeq2 dataset object
dds <- DESeqDataSetFromMatrix(
  countData = count_data,
  colData = sample_info,
  design = ~ Treatment
)

# Attach gene annotations to the DESeq2 object
rowData(dds) <- gene_data

# Perform size factor normalization
dds <- estimateSizeFactors(dds)

# Apply variance stabilizing transformation (VST)
vst_data <- vst(dds, blind = TRUE)
normalized_counts <- assay(vst_data)

# Save normalized expression data
write.table(normalized_counts, "normalized_counts.tsv", sep = "\t", row.names = TRUE, quote = FALSE)





# STEP 5: Identify differentially expressed genes
#-----------------------------------------------

# Perform differential expression analysis
dds <- DESeq(dds)

# Extract results comparing EV_AKO treatment to EV control
res <- results(dds, contrast = c("Treatment", "EV_AKO", "EV"))

# Filter for significant DEGs
# Adjusted p-value < 0.05 and absolute log2 fold change > 1
sig_genes <- res[which(res$padj < 0.05 & abs(res$log2FoldChange) > 1), ]

# Extract DEG names and their expression values
deg_names <- rownames(sig_genes)
deg_expression <- normalized_counts[deg_names, ]



----------------------------------------------
  # STEP 1: Load pre-processed RNA-seq data
  #-----------------------------------------------

# Load normalized expression data
# Rows are genes, columns are samples
expr_data <- read.table("normalized_counts.tsv", header = TRUE, row.names = 1)

# Load sample information
# Contains treatment/condition information for each sample
pheno_data <- read.table("sample_info.tsv", header = TRUE, row.names = 1)

# Load gene annotation data
# Contains ENSEMBL IDs and gene symbols
gene_data <- read.table("gene_data.tsv", header = TRUE)





#-----------------------------------------------
# STEP 2: Quality control and data filtering
#-----------------------------------------------

# Check for missing values
missing_count <- sum(is.na(expr_data))
cat("Number of missing values:", missing_count, "\n")
# Number of missing values: 0

# Remove genes with missing values if any
if (missing_count > 0) {
  expr_data <- expr_data[complete.cases(expr_data), ]
  cat("Removed genes with missing values. New dimension:", dim(expr_data), "\n")
}

# Calculate variance for each gene
gene_variance <- apply(expr_data, 1, var)

# Visualize variance distribution
png("plots/01_gene_variance_distribution.png", width = 800, height = 600, res = 100)
hist(log10(gene_variance), 
     breaks = 50,
     main = "Distribution of Gene Variance (log10)",
     xlab = "log10(Variance)",
     col = "lightblue",
     border = "white")
abline(v = log10(quantile(gene_variance, 0.4)), 
       col = "red", lty = 2, lwd = 2)
legend("topright", 
       legend = "40th percentile cutoff",
       col = "red", lty = 2, lwd = 2,
       bty = "n")
dev.off()

# Remove genes with very low variance (bottom 40%)
# These genes provide little information for co-expression analysis
variance_threshold <- quantile(gene_variance, 0.4)
cat("Variance threshold (40th percentile):", variance_threshold, "\n")
# Variance threshold (40th percentile): 0.04054851

expr_data_filtered <- expr_data[gene_variance > variance_threshold, ]
cat("Genes after variance filtering:", nrow(expr_data_filtered), "\n")
# Genes after variance filtering: 10853

# Transpose for WGCNA (samples as rows, genes as columns)
expr_matrix <- t(expr_data_filtered)

# Check sample clustering to identify outliers
sample_tree <- hclust(dist(expr_matrix), method = "average")
png("plots/02_sample_clustering_outliers.png", width = 1200, height = 600, res = 100)
par(mar = c(2, 5, 2, 2))
plot(sample_tree, 
     main = "Sample Clustering to Detect Outliers",
     sub = "Height represents dissimilarity between samples",
     xlab = "", 
     cex = 0.7)
abline(h = 200, col = "red", lty = 2, lwd = 2)
dev.off()





#-----------------------------------------------
# STEP 3: Prepare trait data matrix
#-----------------------------------------------

# Create binary trait matrix for each experimental condition
# Each column represents presence (1) or absence (0) of that condition
trait_data <- data.frame(
  EV = as.numeric(pheno_data$Treatment == "EV"),
  Y = as.numeric(pheno_data$Treatment == "Y"),
  T = as.numeric(pheno_data$Treatment == "T"),
  N = as.numeric(pheno_data$Treatment == "N"),
  EV_AKO = as.numeric(pheno_data$Treatment == "EV_AKO"),
  Y_AKO = as.numeric(pheno_data$Treatment == "Y_AKO"),
  T_AKO = as.numeric(pheno_data$Treatment == "T_AKO"),
  N_AKO = as.numeric(pheno_data$Treatment == "N_AKO")
)

# Set row names to match expression matrix
rownames(trait_data) <- rownames(expr_matrix)

# Verify that sample names match between expression and trait data
if (!all(rownames(expr_matrix) == rownames(trait_data))) {
  stop("ERROR: Sample names don't match between expression and trait data!")
}


#-----------------------------------------------
# STEP 4: Verify data structure
#-----------------------------------------------

# Visualize expression data distribution
png("plots/03_expression_data_distribution.png", width = 1200, height = 600, res = 100)
par(mfrow = c(1, 2))

# Overall distribution
hist(as.matrix(expr_matrix), 
     breaks = 100, 
     main = "Distribution of Expression Values",
     sub = "Should be roughly normal for VST/rlog data",
     xlab = "Expression Value",
     col = "lightblue",
     border = "white")

# Sample distributions (boxplot)
boxplot(expr_matrix[, sample(1:ncol(expr_matrix), min(20, ncol(expr_matrix)))], 
        main = "Sample Expression Distributions\n(20 random genes)",
        xlab = "Genes",
        ylab = "Expression Value",
        las = 2,
        cex.axis = 0.6,
        col = "lightblue",
        border = "gray40")

par(mfrow = c(1, 1))
dev.off()




###########

disableWGCNAThreads()

sft <- pickSoftThreshold(
  datExpr,
  powerVector = powers,
  verbose = 5
)

#-----------------------------------------------
# STEP 5: Choose soft-thresholding power
#-----------------------------------------------

# Calculate network topology for various soft-thresholding powers
powers <- c(seq(1, 10, by = 1), seq(12, 20, by = 2))

# This function tests different power values
sft <- pickSoftThreshold(
  expr_matrix,
  powerVector = powers,
  networkType = "signed",
  verbose = 3
)

# Plot scale-free topology fit index as a function of power
par(mfrow = c(1, 2))

# Scale-free topology fit
png("plots/04_soft_thresholding_power_selection.png", width = 1200, height = 600, res = 100)
par(mfrow = c(1, 2))

# Scale-free topology fit
plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     main = "Scale Independence",
     type = "n")
text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, 
     cex = 0.9, 
     col = "red")
abline(h = 0.85, col = "red", lty = 2)

# Mean connectivity
plot(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     main = "Mean Connectivity",
     type = "n")
text(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     labels = powers, 
     cex = 0.9, 
     col = "red")

par(mfrow = c(1, 1))
dev.off()

# Choose the power: first value that reaches plateau in scale-free topology
# and gives reasonable mean connectivity
soft_power <- sft$powerEstimate
if (is.na(soft_power)) {
  soft_power <- 6  # Use 6 as default if automatic selection fails
  cat("Warning: Automatic power selection failed. Using default power =", soft_power, "\n")
} else {
  cat("Selected soft-thresholding power:", soft_power, "\n")
}





#-----------------------------------------------
# STEP 6: Construct network and identify modules
#-----------------------------------------------

# CRITICAL: Fix for WGCNA cor() function bug
# WGCNA has its own cor() function that conflicts with stats::cor()
cor <- WGCNA::cor

# One-step network construction and module detection
net <- blockwiseModules(
  expr_matrix,
  power = soft_power,
  networkType = "signed",
  TOMType = "signed",
  minModuleSize = 30,
  reassignThreshold = 0,
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs = TRUE,
  saveTOMFileBase = "GSE261875_TOM",
  verbose = 3
)

# Restore the standard cor() function after network construction
cor <- stats::cor

# Convert numeric module labels to colors for visualization
module_colors <- labels2colors(net$colors)

# Count genes in each module
module_table <- table(module_colors)
cat("\nNumber of modules identified:", length(unique(module_colors)) - 1, "\n")
# Number of modules identified: 27 
print(module_table)
# black blue brown cyan
# 560 1694 1060 216 

# Plot the dendrogram with module colors
png("plots/05_gene_dendrogram_modules.png", width = 1200, height = 800, res = 100)
plotDendroAndColors(
  net$dendrograms[[1]],
  module_colors[net$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05,
  main = "Gene Dendrogram and Module Colors"
)
dev.off()


#-----------------------------------------------
# STEP 7: Calculate and visualize module eigengenes
#-----------------------------------------------

# Extract module eigengenes
module_eigengenes <- net$MEs

# Reorder modules by hierarchical clustering of eigengenes
module_eigengenes <- orderMEs(module_eigengenes)

# Visualize module eigengenes across samples
# Create a heatmap showing how each module's expression varies across samples
png("plots/06_module_eigengenes_heatmap.png", width = 1200, height = 1000, res = 100)
pheatmap(
  t(module_eigengenes),
  scale = "row",
  clustering_distance_cols = "euclidean",
  clustering_method = "average",
  annotation_col = trait_data,
  show_colnames = FALSE,
  main = "Module Eigengenes Across Samples",
  color = colorRampPalette(c("blue", "white", "red"))(50),
  fontsize = 8
)
dev.off()


#-----------------------------------------------
# STEP 8: Module-trait association analysis
#-----------------------------------------------

# Calculate correlations between module eigengenes and traits
module_trait_cor <- cor(module_eigengenes, trait_data, use = "pairwise.complete.obs")

# Calculate p-values for these correlations
nSamples <- nrow(expr_matrix)
module_trait_pvalue <- corPvalueStudent(module_trait_cor, nSamples)

# Create a text matrix for the heatmap
# Format: correlation coefficient (p-value)
textMatrix <- paste(signif(module_trait_cor, 2), "\n(",
                    signif(module_trait_pvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(module_trait_cor)

# Visualize module-trait relationships with a heatmap
png("plots/07_module_trait_relationships.png", width = 1200, height = 1000, res = 100)
par(mar = c(8, 10, 3, 3))
labeledHeatmap(
  Matrix = module_trait_cor,
  xLabels = names(trait_data),
  yLabels = names(module_eigengenes),
  ySymbols = names(module_eigengenes),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  cex.text = 0.5,
  zlim = c(-1, 1),
  main = "Module-Trait Relationships"
)
dev.off()



#-----------------------------------------------
# STEP 9: Gene significance and module membership
#-----------------------------------------------

# 1. Calculate Gene Significance (GS)
# Ensure nSamples is defined (usually nrow(expr_matrix))
nSamples <- nrow(expr_matrix)
gene_trait_cor <- cor(expr_matrix, trait_data$T_AKO, use = "pairwise.complete.obs")
gene_trait_pvalue <- corPvalueStudent(as.numeric(gene_trait_cor), nSamples)

# CRITICAL: Ensure gene names are assigned
names(gene_trait_cor) <- colnames(expr_matrix)
names(gene_trait_pvalue) <- colnames(expr_matrix)

# 2. Identify the strongest module (e.g., ME8)
t_ako_correlations <- abs(module_trait_cor[, "T_AKO"])
t_ako_module <- names(module_eigengenes)[which.max(t_ako_correlations)]
cat("\nModule most strongly associated with T_AKO:", t_ako_module, "\n")

# Convert numeric ME label (e.g., "ME8") to color name (e.g., "pink")
numeric_label <- gsub("ME", "", t_ako_module)
module_name <- labels2colors(as.numeric(numeric_label))
cat("Module color name:", module_name, "\n")

# 3. FIX: Get genes in this module 
# We map the colors to the actual column names of the expression matrix
gene_ids <- colnames(expr_matrix)
module_genes <- gene_ids[module_colors == module_name]

cat("Number of genes in", module_name, "module:", length(module_genes), "\n")

# 4. Calculate Module Membership (MM)
# Use the actual column name from module_eigengenes (t_ako_module = "ME8")
gene_module_membership <- cor(expr_matrix[, module_genes, drop = FALSE],
                              module_eigengenes[, t_ako_module, drop = FALSE],
                              use = "pairwise.complete.obs")

# 5. Create Hub Gene Data Frame
module_gene_info <- data.frame(
  Gene = module_genes,
  ModuleMembership = as.numeric(gene_module_membership),
  GeneSignificance = as.numeric(gene_trait_cor[module_genes]),
  MM_pvalue = corPvalueStudent(as.numeric(gene_module_membership), nSamples),
  GS_pvalue = as.numeric(gene_trait_pvalue[module_genes]),
  stringsAsFactors = FALSE
)

# Sort by Module Membership (identify top Hub Genes)
module_gene_info <- module_gene_info[order(-abs(module_gene_info$ModuleMembership)), ]

# 6. Plotting
# (Your plotting code is mostly fine, but let's ensure text position is dynamic)
png(paste0("plots/08_MM_vs_GS_", module_name, "_module.png"), 
    width = 800, height = 600, res = 100)

plot(abs(module_gene_info$ModuleMembership), 
     abs(module_gene_info$GeneSignificance),
     xlab = paste("Module Membership in", module_name, "module"),
     ylab = "Gene Significance for T_AKO",
     main = paste("MM vs GS in", module_name, "module"),
     pch = 20, 
     col = module_name,
     cex = 1.2)

abline(lm(abs(module_gene_info$GeneSignificance) ~ abs(module_gene_info$ModuleMembership)), 
       col = "red", lwd = 2)

mm_gs_cor <- cor(abs(module_gene_info$ModuleMembership), 
                 abs(module_gene_info$GeneSignificance),
                 use = "pairwise.complete.obs")

# Dynamic text placement
legend("topleft", legend = paste("cor =", round(mm_gs_cor, 3)), bty = "n", cex = 1.2)

dev.off()



#-----------------------------------------------
# STEP 10: Identify hub genes in T_AKO-associated module
#-----------------------------------------------

# Calculate connectivity (intramodular connectivity)
# First, calculate adjacency matrix for genes in this module
adjacency_matrix <- adjacency(
  expr_matrix[, module_genes],
  power = soft_power,
  type = "signed"
)

# Calculate connectivity: sum of connection weights for each gene
# (subtract 1 to exclude self-connection)
connectivity <- rowSums(adjacency_matrix) - 1

# Combine all metrics for hub gene identification
hub_gene_info <- data.frame(
  Gene = module_genes,
  Connectivity = connectivity,
  ModuleMembership = as.numeric(gene_module_membership),
  # Use the vector you created in Step 9 and subset it for module genes
  GeneSignificance = gene_trait_cor[module_genes], 
  stringsAsFactors = FALSE
)

# Sort by connectivity to identify hubs
hub_gene_info <- hub_gene_info[order(-hub_gene_info$Connectivity), ]


# Map ENSEMBL IDs to gene symbols for better interpretation
# Merge with gene annotation data
hub_gene_info_annotated <- merge(hub_gene_info, gene_data, 
                                 by.x = "Gene", by.y = "ENSEMBL",
                                 all.x = TRUE)

# Sort again by connectivity after merge
hub_gene_info_annotated <- hub_gene_info_annotated[order(-hub_gene_info_annotated$Connectivity), ]


# Visualize hub genes
# Plot connectivity vs module membership
png(paste0("plots/09_hub_genes_", module_name, "_module.png"), 
    width = 800, height = 600, res = 100)
plot(hub_gene_info$ModuleMembership, 
     hub_gene_info$Connectivity,
     xlab = "Module Membership",
     ylab = "Connectivity (Intramodular)",
     main = paste("Hub Gene Identification in", module_name, "Module"),
     pch = 20,
     col = ifelse(hub_gene_info$Connectivity > quantile(hub_gene_info$Connectivity, 0.9),
                  "red", "black"))
legend("topright", 
       legend = c("Top 10% connected", "Other genes"),
       col = c("red", "black"),
       pch = 20)
dev.off()




#-----------------------------------------------
# STEP 11: Functional enrichment of module genes
#-----------------------------------------------

# 1. Prepare Gene Symbols
# Ensure we only use non-NA symbols for enrichment
module_gene_symbols <- hub_gene_info_annotated$SYMBOL[!is.na(hub_gene_info_annotated$SYMBOL)]

# 2. Perform GO enrichment (Biological Process)
go_enrichment <- enrichGO(
  gene          = module_gene_symbols,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

# Plot GO Enrichment Barplot
if (!is.null(go_enrichment) && nrow(go_enrichment) > 0) {
  png(paste0("plots/11_GO_enrichment_barplot_", module_name, ".png"), 
      width = 1000, height = 800, res = 100)
  # Using print() ensures the plot is sent to the png device
  print(barplot(go_enrichment, 
                showCategory = 15, 
                title = paste("GO Enrichment in", module_name, "Module")))
  dev.off()
}

# 3. KEGG Pathway Enrichment
# First: Convert symbols to Entrez IDs (required for KEGG)
gene_entrez <- bitr(module_gene_symbols, 
                    fromType = "SYMBOL", 
                    toType   = "ENTREZID", 
                    OrgDb    = org.Hs.eg.db)

# Second: Run the enrichment
kegg_enrichment <- enrichKEGG(
  gene          = gene_entrez$ENTREZID,
  organism      = "hsa", # 'hsa' for Homo sapiens
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

# Third: Make results readable (convert Entrez back to Symbols in the table)
if (!is.null(kegg_enrichment) && nrow(kegg_enrichment) > 0) {
  kegg_enrichment <- setReadable(kegg_enrichment, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  
  # Plot KEGG Enrichment Dotplot
  png(paste0("plots/12_KEGG_enrichment_", module_name, ".png"), 
      width = 1000, height = 800, res = 100)
  print(dotplot(kegg_enrichment, 
                showCategory = 15, 
                title = paste("KEGG Pathway Enrichment in", module_name, "Module")))
  dev.off()
} else {
  cat("No significant KEGG pathways found for module:", module_name, "\n")
}



#-----------------------------------------------
# STEP 12: Construct PPI network for hub genes
#-----------------------------------------------
library(STRINGdb)
library(igraph)

# 1. Get top 20 hub genes
# Ensure we handle potential NAs and use the sorted list from earlier
top_hub_genes <- head(hub_gene_info_annotated$SYMBOL[!is.na(hub_gene_info_annotated$SYMBOL)], 20)

# 2. Initialize STRING database 
# Note: score_threshold 400 is "Medium Confidence"
string_db <- STRINGdb$new(version = "11.5", 
                          species = 9606, 
                          score_threshold = 400,
                          input_directory = "") # Leaves directory default if not specified

# 3. Map gene symbols to STRING IDs
hub_genes_df <- data.frame(gene = top_hub_genes, stringsAsFactors = FALSE)
hub_genes_mapped <- string_db$map(hub_genes_df, "gene", removeUnmappedRows = TRUE)

if (nrow(hub_genes_mapped) > 1) {
  # 4. Get interactions
  interactions <- string_db$get_interactions(hub_genes_mapped$STRING_id)
  
  if (nrow(interactions) > 0) {
    cat("Found", nrow(interactions), "interactions among hub genes\n")
    
    # Create mapping for labels
    id_to_gene <- setNames(hub_genes_mapped$gene, hub_genes_mapped$STRING_id)
    
    # Convert STRING IDs to Gene Symbols in the interaction table
    interactions$from_gene <- id_to_gene[interactions$from]
    interactions$to_gene   <- id_to_gene[interactions$to]
    
    # 5. Create igraph object
    # Subset to only unique interactions to avoid double edges in the plot
    ppi_graph <- graph_from_data_frame(
      unique(interactions[, c("from_gene", "to_gene")]),
      directed = FALSE
    )
    
    # 6. Plotting
    png(paste0("plots/13_PPI_network_", module_name, ".png"), 
        width = 1000, height = 1000, res = 120)
    
    par(mar = c(1, 1, 4, 1))
    
    # Use degree to scale node size (more connected = bigger node)
    v_size <- degree(ppi_graph)
    v_size_scaled <- (v_size / max(v_size) * 15) + 10 # Scale between 10-25
    
    plot(ppi_graph,
         vertex.size = v_size_scaled,
         vertex.color = "gold", # Highlighting hubs in gold
         vertex.label.cex = 1.0,
         vertex.label.font = 2, # Bold text
         vertex.label.color = "black",
         edge.color = "gray80",
         edge.width = 1.5,
         layout = layout_with_kk(ppi_graph), # Kamada-Kawai layout often looks better for small PPIs
         main = paste("PPI Network: Top Hubs in", module_name, "Module"))
    
    dev.off()
    
    # 7. Reporting
    cat("\nPPI Network properties:\n")
    cat("  Nodes:", vcount(ppi_graph), "\n")
    cat("  Edges:", ecount(ppi_graph), "\n")
    
    ppi_degree <- degree(ppi_graph)
    cat("\nMost connected proteins (Hubs) in this sub-network:\n")
    print(sort(ppi_degree, decreasing = TRUE))
    
  } else {
    cat("No internal interactions found at this confidence threshold (400).\n")
  }
} else {
  cat("Insufficient genes mapped to STRING database.\n")
}


#-----------------------------------------------
# STEP 13: Global Network Visualization
#-----------------------------------------------
library(WGCNA)
library(igraph)
library(scales) # for alpha transparency

# 1. Select Genes for Visualization
# Using the top 1000 most variable genes to keep the network readable and computational load low
top_genes_var <- names(sort(apply(expr_matrix, 2, var), decreasing = TRUE)[1:1000])

cat("Calculating TOM for top 1000 variable genes...\n")

# 2. Calculate Adjacency and TOM
adj_subset <- adjacency(expr_matrix[, top_genes_var], 
                        power = soft_power, 
                        type = "signed")

tom_subset <- TOMsimilarity(adj_subset)
colnames(tom_subset) <- top_genes_var
rownames(tom_subset) <- top_genes_var

# 3. Apply Thresholding
# We only want to visualize the strongest 5% of connections
tom_threshold <- quantile(tom_subset[upper.tri(tom_subset)], 0.95)
cat("TOM threshold (95th percentile):", tom_threshold, "\n")

adjacency_threshold <- tom_subset
adjacency_threshold[adjacency_threshold < tom_threshold] <- 0

# 4. Create igraph Object
network_graph <- graph_from_adjacency_matrix(
  adjacency_threshold,
  mode = "undirected",
  weighted = TRUE,
  diag = FALSE
)

# 5. Handle Node Metadata (Colors and Sizes)
# IMPORTANT: Ensure module_colors has gene names for correct mapping
names(module_colors) <- colnames(expr_matrix)
node_colors <- module_colors[match(top_genes_var, names(module_colors))]

# Remove isolated nodes (genes with 0 connections after thresholding)
# This significantly cleans up the final visualization
isolated_nodes <- V(network_graph)[degree(network_graph) == 0]
network_graph <- delete.vertices(network_graph, isolated_nodes)
node_colors <- node_colors[-as.numeric(isolated_nodes)]

# Calculate node sizes based on connectivity (hubs appear larger)
node_sizes <- log10(degree(network_graph) + 1) * 4 + 2

# 6. Plotting
png("plots/14_coexpression_network.png", width = 1200, height = 1200, res = 150)
par(mar = c(1, 1, 3, 1), bg = "white")

set.seed(123) # Ensure layout is reproducible
plot(network_graph,
     vertex.size = node_sizes,
     vertex.label = NA,
     vertex.color = node_colors,
     vertex.frame.color = "white",
     vertex.frame.width = 0.2,
     edge.width = 0.4,
     edge.color = alpha("gray80", 0.4),
     layout = layout_with_fr(network_graph), # Fruchterman-Reingold layout
     main = "WGCNA Co-expression Network\n(Top 1000 Variable Genes, Top 5% TOM Connections)")

# 7. Add Legend for the Top Modules
# Identifies the unique colors in this specific subset (excluding grey)
active_colors <- unique(node_colors[node_colors != "grey"])
# Display up to the top 15 modules found in this subset
legend_count <- min(length(active_colors), 15)
legend("topright",
       legend = active_colors[1:legend_count],
       col = active_colors[1:legend_count],
       pch = 19,
       pt.cex = 1.5,
       cex = 0.8,
       title = "Modules",
       bty = "n")

dev.off()

cat("Network visualization saved to plots/14_coexpression_network.png\n")




#-----------------------------------------------
# STEP 14: Module eigengene network
#-----------------------------------------------

# Calculate correlation between module eigengenes
me_correlation <- cor(module_eigengenes, use = "pairwise.complete.obs")

# Calculate p-values
me_pvalue <- corPvalueStudent(me_correlation, nSamples)

# Create adjacency matrix: keep correlations > 0.5 and p < 0.05
me_adjacency <- (abs(me_correlation) > 0.5) & (me_pvalue < 0.05)
diag(me_adjacency) <- FALSE  # Remove self-connections

# Check if there are any edges
n_edges <- sum(me_adjacency) / 2
cat("\nNumber of module pairs with |correlation| > 0.5 and p < 0.05:", n_edges, "\n")
# Number of module pairs with |correlation| > 0.5 and p < 0.05: 106 

if (n_edges == 0) {
  cat("\nNo strong correlations between modules at |r| > 0.5.\n")
  cat("Trying |r| > 0.3...\n\n")
  me_adjacency <- (abs(me_correlation) > 0.3) & (me_pvalue < 0.05)
  diag(me_adjacency) <- FALSE
  n_edges <- sum(me_adjacency) / 2
  correlation_threshold <- 0.3
} else {
  correlation_threshold <- 0.5
}

if (n_edges == 0) {
  cat("No significant correlations between modules even at |r| > 0.3.\n")
  cat("Modules are largely independent.\n")
} else {
  # Create network of module relationships
  me_network <- graph_from_adjacency_matrix(
    me_adjacency,
    mode = "undirected"
  )
  
  # Convert numeric module labels to color names
  module_names_numeric <- gsub("ME", "", names(module_eigengenes))
  module_names_colors <- labels2colors(as.numeric(module_names_numeric))
  
  cat("Module numeric labels:", head(module_names_numeric), "\n")
  cat("Corresponding colors:", head(module_names_colors), "\n")
  
  # Set vertex colors using the actual color names
  V(me_network)$color <- module_names_colors
  
  # Calculate node sizes based on module size
  # Use color names to get counts from module_colors table
  module_sizes_table <- table(module_colors)
  module_sizes_net <- module_sizes_table[module_names_colors]
  
  # Handle any missing values (modules with no genes)
  module_sizes_net[is.na(module_sizes_net)] <- 10
  
  node_sizes_me <- log10(as.numeric(module_sizes_net) + 1) * 10 + 10
  
  # Plot module network
  if (n_edges > 0) {
    png("plots/15_module_eigengene_network.png", width = 1000, height = 1000, res = 100)
    par(mar = c(1, 1, 3, 1))
    set.seed(123)
    plot(me_network,
         vertex.size = node_sizes_me,
         vertex.label = module_names_colors,
         vertex.label.cex = 0.7,
         vertex.label.color = "black",
         vertex.color = module_names_colors,
         vertex.frame.color = "black",
         edge.width = 3,
         edge.color = "gray40",
         layout = layout_with_fr(me_network),
         main = paste0("Module Eigengene Network\n(|correlation| > ", 
                       correlation_threshold, ", p < 0.05)"))
    dev.off()
  }}


# Also show a correlation heatmap for reference
# Use color names for row/column labels
me_cor_labeled <- me_correlation
rownames(me_cor_labeled) <- paste0(module_names_colors, " (", gsub("ME", "", rownames(me_cor_labeled)), ")")
colnames(me_cor_labeled) <- paste0(module_names_colors, " (", gsub("ME", "", colnames(me_cor_labeled)), ")")

png("plots/16_module_eigengene_correlation_heatmap.png", width = 1200, height = 1200, res = 100)
pheatmap(
  me_cor_labeled,
  color = colorRampPalette(c("blue", "white", "red"))(50),
  breaks = seq(-1, 1, length.out = 51),
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  main = "Module Eigengene Correlation Heatmap",
  fontsize = 7
)
dev.off()




#-----------------------------------------------
# STEP 15: Module expression heatmap
#-----------------------------------------------

# Select the T_AKO-associated module
module_numeric <- gsub("ME", "", t_ako_module)
cat("Module numeric label:", module_numeric, "\n")

# Convert numeric label to color name
module_of_interest <- labels2colors(as.numeric(module_numeric))
cat("Module color name:", module_of_interest, "\n")

# Get indices of genes in this module
module_gene_indices <- which(module_colors == module_of_interest)
cat("Number of genes in", module_of_interest, "module:", length(module_gene_indices), "\n")

# Check if we found genes
if (length(module_gene_indices) == 0) {
  cat("\nERROR: No genes found in module", module_of_interest, "\n")
  cat("Available colors in module_colors:", unique(module_colors), "\n")
  stop("Cannot create heatmap - no genes in selected module")
}

# If module is too large, select top genes by module membership
max_genes_for_heatmap <- 100
if (length(module_gene_indices) > max_genes_for_heatmap) {
  # Get MM values for genes in this module
  mm_values <- abs(cor(expr_matrix[, module_gene_indices],
                       module_eigengenes[, t_ako_module],
                       use = "pairwise.complete.obs"))
  
  # Select top genes by MM
  top_mm_indices <- order(mm_values, decreasing = TRUE)[1:max_genes_for_heatmap]
  selected_genes <- colnames(expr_matrix)[module_gene_indices[top_mm_indices]]
  cat("Selecting top", max_genes_for_heatmap, "genes by module membership\n")
} else {
  selected_genes <- colnames(expr_matrix)[module_gene_indices]
}

# Create annotation for samples
sample_annotation <- trait_data[, c("T_AKO", "EV", "EV_AKO")]
colnames(sample_annotation) <- c("T_AKO", "Control", "AKO_Background")

# Create heatmap
png(paste0("plots/17_expression_heatmap_", module_of_interest, "_module.png"), 
    width = 1200, height = 1000, res = 100)
pheatmap(
  t(expr_matrix[, selected_genes]),
  scale = "row",
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "euclidean",
  annotation_col = sample_annotation,
  show_rownames = FALSE,
  show_colnames = FALSE,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  main = paste("Expression Heatmap:", module_of_interest, "Module\n",
               "(Top", length(selected_genes), "genes by Module Membership)"),
  annotation_colors = list(
    T_AKO = c("0" = "white", "1" = "darkred"),
    Control = c("0" = "white", "1" = "darkgreen"),
    AKO_Background = c("0" = "white", "1" = "darkblue")
  )
)
dev.off()

# Save module gene list with annotations
module_gene_list <- data.frame(
  ENSEMBL = selected_genes,
  Module = module_of_interest,
  ModuleMembership = as.numeric(cor(expr_matrix[, selected_genes],
                                    module_eigengenes[, t_ako_module],
                                    use = "pairwise.complete.obs")),
  stringsAsFactors = FALSE
)

# Add gene symbols
module_gene_list <- merge(module_gene_list, gene_data, 
                          by = "ENSEMBL", all.x = TRUE)
module_gene_list <- module_gene_list[order(-abs(module_gene_list$ModuleMembership)), ]

# Save to file
output_file <- paste0(module_of_interest, "_module_genes.txt")
write.table(module_gene_list,
            file = output_file,
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)






#-----------------------------------------------
# STEP 16: Build networks for WT and AKO conditions
#-----------------------------------------------

# Separate samples by genetic background
wt_samples <- rownames(trait_data)[trait_data$EV == 1 | 
                                     trait_data$Y == 1 | 
                                     trait_data$T == 1 | 
                                     trait_data$N == 1]
ako_samples <- rownames(trait_data)[trait_data$EV_AKO == 1 | 
                                      trait_data$Y_AKO == 1 | 
                                      trait_data$T_AKO == 1 | 
                                      trait_data$N_AKO == 1]


# Subset expression matrices
expr_wt <- expr_matrix[wt_samples, ]
expr_ako <- expr_matrix[ako_samples, ]

# Build network for WT condition
cor <- WGCNA::cor
net_wt <- blockwiseModules(
  expr_wt, 
  power = soft_power,
  networkType = "signed",
  TOMType = "signed",
  minModuleSize = 30,
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  verbose = 0
)
cor <- stats::cor

# Build network for AKO condition
cor <- WGCNA::cor
net_ako <- blockwiseModules(
  expr_ako, 
  power = soft_power,
  networkType = "signed",
  TOMType = "signed",
  minModuleSize = 30,
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  verbose = 0
)
cor <- stats::cor

# Convert to color labels
colors_wt <- labels2colors(net_wt$colors)
colors_ako <- labels2colors(net_ako$colors)

cat("\nWT network: ", length(unique(colors_wt)) - 1, "modules\n")
# WT network:  17 modules
cat("AKO network:", length(unique(colors_ako)) - 1, "modules\n")
# AKO network: 12 modules

# Plot 18a: WT Dendrogram
png("plots/18a_dendrogram_WT.png", width = 1200, height = 600, res = 100)
par(mar = c(2, 5, 3, 2))
plotDendroAndColors(
  net_wt$dendrograms[[1]],
  colors_wt[net_wt$blockGenes[[1]]],
  "WT Modules",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05,
  main = "Gene Dendrogram - Wild-Type Condition"
)
dev.off()

# Plot 18b: AKO Dendrogram
png("plots/18b_dendrogram_AKO.png", width = 1200, height = 600, res = 100)
par(mar = c(2, 5, 3, 2))
plotDendroAndColors(
  net_ako$dendrograms[[1]],
  colors_ako[net_ako$blockGenes[[1]]],
  "AKO Modules",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05,
  main = "Gene Dendrogram - AKO Condition"
)
dev.off()



#-----------------------------------------------
# STEP 17: Module preservation analysis
#-----------------------------------------------

# CRITICAL: Check for genes with zero variance in each subset
# This is required before modulePreservation

# Check WT samples
gsg_wt <- goodSamplesGenes(expr_wt, verbose = 0)
cat("WT - Good samples:", sum(gsg_wt$goodSamples), "/", length(gsg_wt$goodSamples), "\n")
# WT - Good samples: 16 / 16
cat("WT - Good genes:", sum(gsg_wt$goodGenes), "/", length(gsg_wt$goodGenes), "\n")
# WT - Good genes: 10847 / 10853 

# Check AKO samples
gsg_ako <- goodSamplesGenes(expr_ako, verbose = 0)
cat("AKO - Good samples:", sum(gsg_ako$goodSamples), "/", length(gsg_ako$goodSamples), "\n")
# AKO - Good samples: 16 / 16 
cat("AKO - Good genes:", sum(gsg_ako$goodGenes), "/", length(gsg_ako$goodGenes), "\n")
# AKO - Good genes: 10853 / 10853 

# Keep only genes that are good in BOTH conditions
genes_to_keep <- gsg_wt$goodGenes & gsg_ako$goodGenes
cat("\nGenes passing QC in both conditions:", sum(genes_to_keep), "/", length(genes_to_keep), "\n")

if (sum(!genes_to_keep) > 0) {
  cat("Removing", sum(!genes_to_keep), "genes with zero variance or missing data\n\n")
  # Genes passing QC in both conditions: 10847 / 10853 
  
  # Filter expression matrices
  expr_wt_filtered <- expr_wt[, genes_to_keep]
  expr_ako_filtered <- expr_ako[, genes_to_keep]
  
  # Filter module colors to match
  colors_wt_filtered <- colors_wt[genes_to_keep]
  colors_ako_filtered <- colors_ako[genes_to_keep]
  
} else {
  expr_wt_filtered <- expr_wt
  expr_ako_filtered <- expr_ako
  colors_wt_filtered <- colors_wt
  colors_ako_filtered <- colors_ako
}

# Verify dimensions match
cat("Filtered WT dimensions:", dim(expr_wt_filtered), "\n")
# Filtered WT dimensions: 16 10847 
cat("Filtered AKO dimensions:", dim(expr_ako_filtered), "\n")
# Filtered AKO dimensions: 16 10847 

# Prepare data for preservation analysis
multiExpr <- list(
  WT = list(data = expr_wt_filtered),
  AKO = list(data = expr_ako_filtered)
)

multiColor <- list(
  WT = colors_wt_filtered
)

# Calculate preservation statistics
# Compare WT modules in AKO network

set.seed(123)  # For reproducibility
preservation <- modulePreservation(
  multiExpr,
  multiColor,
  referenceNetworks = 1,  # WT is reference
  nPermutations = 100,     # Use 200+ for publication-quality results
  verbose = 3
)

# Extract preservation statistics
preservation_stats <- preservation$preservation$Z$ref.WT$inColumnsAlsoPresentIn.AKO

# Add module names
preservation_stats$module <- rownames(preservation_stats)
preservation_stats$color <- rownames(preservation_stats)

# Classify preservation status
# Zsummary > 10: strong preservation
# Zsummary 2-10: moderate preservation
# Zsummary < 2: no preservation
preservation_stats$status <- cut(
  preservation_stats$Zsummary.pres,
  breaks = c(-Inf, 2, 10, Inf),
  labels = c("Not Preserved", "Moderate", "Strong")
)


# Visualize preservation
# Plot 19: Module Preservation
png("plots/19_module_preservation.png", width = 1200, height = 800, res = 100)
ggplot(preservation_stats, 
       aes(x = moduleSize, 
           y = Zsummary.pres,
           color = status,
           label = color)) +
  geom_point(size = 4, alpha = 0.7) +
  geom_hline(yintercept = c(2, 10), 
             linetype = "dashed", 
             color = "gray50",
             linewidth = 0.5) +
  geom_text_repel(size = 3, max.overlaps = 20, 
                  segment.color = "gray70",
                  box.padding = 0.5) +
  scale_color_manual(values = c("Not Preserved" = "#d62728",
                                "Moderate" = "#ff7f0e",
                                "Strong" = "#2ca02c")) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11, color = "gray40")
  ) +
  labs(
    title = "Module Preservation Analysis",
    subtitle = "WT modules tested in AKO network",
    x = "Module Size (number of genes)",
    y = "Preservation Z-summary",
    color = "Preservation Status",
    caption = "Zsummary > 10: Strong | 2-10: Moderate | <2: Not preserved"
  )
dev.off()

# Save preservation results
write.table(preservation_stats,
            file = "module_preservation_results.txt",
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)





#-----------------------------------------------
# STEP 18: Compare module eigengenes across conditions
#-----------------------------------------------

# We'll compare the T_AKO-associated module between WT and AKO backgrounds
# Using the module assignments from our combined analysis (Step 6)

# Convert numeric module label to color name FIRST
module_numeric <- gsub("ME", "", t_ako_module)
module_color <- labels2colors(as.numeric(module_numeric))
module_to_compare <- paste0("ME", module_color)

# Calculate module eigengenes for WT samples
# Using the same module assignments from the combined analysis
me_wt <- moduleEigengenes(expr_wt, 
                          colors = module_colors[colnames(expr_wt)])$eigengenes

# Calculate module eigengenes for AKO samples  
# Using the same module assignments from the combined analysis
me_ako <- moduleEigengenes(expr_ako, 
                           colors = module_colors[colnames(expr_ako)])$eigengenes


# Create comparison data frame
me_comparison <- data.frame(
  Eigengene = c(me_wt[, module_to_compare], me_ako[, module_to_compare]),
  Background = c(rep("WT", nrow(me_wt)), rep("AKO", nrow(me_ako))),
  Treatment = c(pheno_data[wt_samples, "Treatment"],
                gsub("_AKO", "", pheno_data[ako_samples, "Treatment"]))
)


# Plot 20: Overall comparison (WT vs AKO)
png(paste0("plots/20_ME_comparison_overall_", module_color, ".png"), 
    width = 800, height = 600, res = 100)
ggplot(me_comparison, aes(x = Background, y = Eigengene, fill = Background)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 2.5) +
  theme_minimal(base_size = 12) +
  scale_fill_manual(values = c("WT" = "#3498db", "AKO" = "#e74c3c")) +
  labs(
    title = paste("Module Eigengene Comparison:", module_color, "Module"),
    subtitle = "Overall expression pattern between backgrounds",
    y = "Module Eigengene Value",
    x = "Genetic Background"
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10, color = "gray40")
  )
dev.off()

# Plot 21: Detailed comparison by treatment
png(paste0("plots/21_ME_comparison_by_treatment_", module_color, ".png"), 
    width = 1000, height = 600, res = 100)
ggplot(me_comparison, aes(x = Treatment, y = Eigengene, fill = Background)) +
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7) +
  geom_point(position = position_dodge(width = 0.8), alpha = 0.6, size = 2) +
  theme_minimal(base_size = 12) +
  scale_fill_manual(values = c("WT" = "#3498db", "AKO" = "#e74c3c")) +
  labs(
    title = paste("Module Eigengene by Treatment:", module_color, "Module"),
    subtitle = "Expression across treatments and backgrounds",
    y = "Module Eigengene Value",
    x = "Treatment",
    fill = "Background"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10, color = "gray40"),
    legend.position = "bottom"
  )
dev.off()
