library(RWeka)
library(mice)
library(Hmisc)
library(VIM)
library(readxl)
library(ggplot2)
library(openxlsx)
library(dplyr)
library(UpSetR)
library(writexl)
library(pheatmap)
library(glmnet)

setwd("/Users/mortezaabyadeh/Desktop")
dataset <- read_excel("set5.xlsx", sheet = "DAPs")
head(dataset)

md.pattern(dataset)

mice_plot <- aggr(dataset, col=c("green","red"),
                  numbers=TRUE, sortVars=TRUE,
                  labels=names(dataset), cex.axis=.7,
                  gap=3, ylab=c("Missing data","Pattern"))
dim(dataset)

mrnaNorm <- read.table("dataset.txt", 
                       header = F, fill = T, skip = 1)
mrnaIDs <- read.table("dataset.txt", 
                      header = F, fill = T, nrows = 1)

mrnaIDs <- mrnaIDs[,-1]
colnames(mrnaIDs)

samp <- lapply(as.list(t(mrnaIDs)), function(t) substr(unlist(strsplit(t, "_"))[1], 1, 3))
sampleType <- as.data.frame(samp)
head(sampleType)
colnames(sampleType)

sampClass <- lapply(samp, function(t) if (t == "Eto") return(1) else return(0))
mrnaClass <- as.data.frame(sampClass)
head(mrnaClass)

table(unlist(sampClass))

sampClassNum <- lapply(samp, function(t) if (t == "Eto") return(1) else return(0))
mrnaClassNum <- as.data.frame(sampClassNum)

head(mrnaNorm, 3)

geneNames <- mrnaNorm[1]  
dim(geneNames)
head(geneNames)

colnames(mrnaNorm)[1:3]
mrnaData = t(mrnaNorm[, -1])

# Prepare data for Lasso
X <- as.matrix(mrnaData)
y <- as.numeric(unlist(mrnaClassNum))

# Apply Lasso for feature selection ########################## I used k fold due to low nuber of samples #####
cv_lasso <- cv.glmnet(X, y, alpha = 1, family = "binomial")

# Get the indices of the most relevant genes
lasso_coef <- coef(cv_lasso, s = "lambda.min")
selected_genes_indices <- which(lasso_coef != 0)[-1]  # Exclude the intercept

selected_gene_names <- geneNames[, 1][selected_genes_indices]



####################################################### used k fold starting from here, do not use the above one


set.seed(42)
cv_lasso <- cv.glmnet(X, y, alpha = 1, family = "binomial", nfolds = 6)


lasso_coefs_matrix <- as.matrix(lasso_coefs)

# Extract the gene names corresponding to non-zero coefficients
selected_genes <- rownames(lasso_coefs_matrix)[lasso_coefs_matrix != 0]

# Display selected genes
print(selected_genes)

plot(cv_lasso)
title("Cross-Validation Results for Lasso")
abline(v = log(cv_lasso$lambda.min), col = "red", lty = 2)  # for lambda.min
abline(v = log(cv_lasso$lambda.1se), col = "blue", lty = 2) 

############################################### again same results due to small sample size #############


mrnaDataReduced <- mrnaData[, selected_genes_indices]
length(selected_gene_names)

transpose_mrnaDataReduced <- t(mrnaDataReduced)
dim(transpose_mrnaDataReduced)

colnames(transpose_mrnaDataReduced) <- c(rep("Eto", 3), rep("FTY_Eto", 3))

selected_genes_expression <- cbind(Gene = selected_gene_names, transpose_mrnaDataReduced)
expression_df <- as.data.frame(selected_genes_expression)

head(expression_df)

write_xlsx(expression_df, "selected_genes_expression_lasso.xlsx")

proteome <- read_excel("set5.xlsx", sheet="DAPs")
proteome <- proteome[, c("Gene", "p-adj", "log2FC")]

dataset1_2 <- merge(proteome, expression_df, by = "Gene")


upset(fromList(list(Set1 = dataset1_2$Gene, Set2 = expression_df$Gene)))

write_xlsx(dataset1_2, "selected_genes_expression_final.xlsx")
