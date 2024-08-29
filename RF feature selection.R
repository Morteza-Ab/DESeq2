library(randomForest)
library(readxl)
library(writexl)
library(dplyr)
library(UpSetR)

setwd("/Users/mortezaabyadeh/Desktop")
dataset <- read_excel("set5.xlsx", sheet = "DAPs")
mrnaNorm <- read.table("dataset.txt", header = F, fill = T, skip = 1)
mrnaIDs <- read.table("dataset.txt", header = F, fill = T, nrows = 1)
print(mrnaNorm)

mrnaIDs <- mrnaIDs[, -1]
colnames(mrnaIDs)
samp <- lapply(as.list(t(mrnaIDs)), function(t) substr(unlist(strsplit(t, "_"))[1], 1, 3))
sampClassNum <- lapply(samp, function(t) if (t == "Eto") return(1) else return(0))
mrnaClassNum <- as.data.frame(sampClassNum)
geneNames <- mrnaNorm[1]  
mrnaData <- t(mrnaNorm[, -1])

X <- as.matrix(mrnaData)
y <- as.factor(unlist(mrnaClassNum))

rf_model <- randomForest(X, y, importance = TRUE, ntree = 100)

# Extract feature importance
importance_scores <- importance(rf_model)

# Order features by importance
top_scores <- order(importance_scores[, "MeanDecreaseAccuracy"], decreasing = TRUE)

top_n <- 30
top_feature_indices <- top_scores[1:top_n]

top_feature_names <- rownames(importance_scores)[top_feature_indices]
top_feature_importance <- importance_scores[top_feature_indices, ]

print(top_feature_names)
print(top_feature_importance)


selected_gene_names <- geneNames[, 1][top_features_indices]

mrnaDataReduced <- mrnaData[, top_features_indices]
transpose_mrnaDataReduced <- t(mrnaDataReduced)
dim(transpose_mrnaDataReduced)
print(selected_gene_names)
length(selected_gene_names)
dim(transpose_mrnaDataReduced)

colnames(transpose_mrnaDataReduced) <- c(rep("Eto", 3), rep("FTY_Eto",3))
head(transpose_mrnaDataReduced)

transpose_mrnaDataReduced <- cbind(top_feature_importance, transpose_mrnaDataReduced)
print(transpose_mrnaDataReduced)
selected_genes_expression <- cbind(Gene = selected_gene_names, transpose_mrnaDataReduced)
expression_df <- as.data.frame(selected_genes_expression)
head(expression_df)
write_xlsx(expression_df, "selected_genes_expression_rf.xlsx")

proteome <- read_excel("set5.xlsx", sheet = "DAPs")
proteome <- proteome[, c("Gene", "p-adj", "log2FC")]
dataset1_2 <- merge(proteome, expression_df, by = "Gene")

upset(fromList(list(Set1 = dataset1_2$Gene, Set2 = expression_df$Gene)))

write_xlsx(dataset1_2, "selected_genes_expression_final_rf.xlsx")

barplot(importance_scores[top_features_indices, "MeanDecreaseAccuracy"], 
        names.arg = selected_gene_names, 
        las = 2, 
        main = "Top 30 Features Importance")
