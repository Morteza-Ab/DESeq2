library(shiny)
library(ggplot2)
library(reshape2)
install.packages("ggcorrplot")
library(grid)
library(gridExtra)
library(ggcorrplot)
library(DESeq2)
library(readxl)
library(dplyr)
library(pheatmap)
library(data.table)
library(ggplot2)
library(plotly)
library(htmlwidgets)
library(RColorBrewer) 
library(ggrepel)
library(tidyverse)
library(reshape2)
library(DOSE)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(enrichR)
library(biomaRt)
library(RWeka)
library(mice)
library(Hmisc)
library(VIM)
library(UpSetR)
library(writexl)

transcriptome_Eto_vs_CTRL <- read.xlsx("https://github.com/Morteza-Ab/Shiny-Senescence/raw/main/clean_transcriptome_Eto_vs_CTRL.xlsx")
transcriptome_FTY_vs_CTRL <- read.xlsx("https://github.com/Morteza-Ab/Shiny-Senescence/raw/main/clean_transcriptome_FTY_vs_CTRL.xlsx")
transcriptome_FTY_Eto_vs_CTRL <- read.xlsx("https://github.com/Morteza-Ab/Shiny-Senescence/raw/main/clean_transcriptome_FTY_Eto_vs_CTRL.xlsx")
transcriptome_FTY_vs_Eto <- read.xlsx("https://github.com/Morteza-Ab/Shiny-Senescence/raw/main/clean_transcriptome_FTY_vs_Eto.xlsx")
transcriptome_FTY_Eto_vs_Eto <- read.xlsx("https://github.com/Morteza-Ab/Shiny-Senescence/raw/main/clean_transcriptome_FTY_Eto_vs_Eto.xlsx")
transcriptome_FTY_Eto_vs_FTY <- read.xlsx("https://github.com/Morteza-Ab/Shiny-Senescence/raw/main/clean_transcriptome_FTY_Eto_vs_FTY.xlsx")


head(transcriptome_Eto_vs_CTRL)
class(transcriptome_Eto_vs_CTRL)
head(transcriptome_FTY_vs_Eto)
head(transcriptome_FTY_vs_CTRL)
head(transcriptome_FTY_Eto_vs_CTRL)
head(transcriptome_FTY_Eto_vs_Eto)
head(transcriptome_FTY_Eto_vs_FTY)


# List of all datasets
datasets <- list(
  transcriptome_Eto_vs_CTRL,
  transcriptome_FTY_vs_Eto,
  transcriptome_FTY_vs_CTRL,
  transcriptome_FTY_Eto_vs_CTRL,
  transcriptome_FTY_Eto_vs_Eto,
  transcriptome_FTY_Eto_vs_FTY
)

dim(transcriptome_FTY_vs_CTRL)
dim(transcriptome_FTY_Eto_vs_FTY)

# Extract the common genes
common_genes <- Reduce(intersect, lapply(datasets, function(df) df$Gene))

# Subset each dataset to include only the common genes and first 7 columns
datasets_subset <- lapply(datasets, function(df) {
  df <- df[df$Gene %in% common_genes, ]
  df[, 1:7] # Keep only the first 7 columns
})

# Merge datasets by "Gene"
merged_data <- Reduce(function(x, y) merge(x, y, by = "Gene"), datasets_subset)

# View the result
head(merged_data)

dim(merged_data)

data <- merged_data[,c("Gene", "CTRL-1", "CTRL-2", "CTRL-3", "Eto-1",  "Eto-2",  "Eto-3", "FTY-1", "FTY-2",  "FTY-3","FTY-Eto-1", "FTY-Eto-2", "FTY-Eto-3")]

colnames(transcriptome_FTY_Eto_vs_Eto)
merged_data[2:4, ]

dim(data)
head(data)


head(data1)


# Install and load the openxlsx package
if (!require("openxlsx")) install.packages("openxlsx")
library(openxlsx)
getwd()
# Write data1 to an Excel file
write.xlsx(data1, file = "/Users/mortezaabyadeh/Documents/npc-fty/whole_RNAseq.xlsx", rowNames = FALSE)












data[,1] <- rownames(data)
data <- data[,-1]

data <- na.omit(data)
data <- data.frame(lapply(data, function(x) if(is.numeric(x)) as.integer(x) else x))
head(data)
gr <- c(rep("CTRL", 3), rep("Eto", 3), rep("FTY", 3), rep("FTY-Eto", 3))
gr <- as.factor(gr)

dim(data)
head(data, 5)
data <- na.omit(data)
head(data)

colnames(data) <- gr
colData <- data.frame(Group = gr)
head(data)

rownames(colData) <- colnames(data)

cds <- DESeqDataSetFromMatrix(countData = data, colData = colData, design = ~ Group)
cds
cds <- DESeq(cds)


data.norm <- log2(1+counts(cds, normalized=T))
head(data.norm)
boxplot(data.norm)

data.mean.center <- t(scale(t(data.norm), scale = F))

### PCA

pc <- prcomp(data.mean.center)
pc1 <- prcomp(data.norm)
plot(pc1)
plot(pc)
head(pc$rotation)
pcr <- data.frame(pc$r)
pcr1 <- data.frame(pc1$r)
head(pcr)
head(data)
pcr1$group <- gr
pcr$group <- gr
head(pcr)


pc_var <- pc$sdev^2 / sum(pc$sdev^2)
pc_var

pc1_var <- pc1$sdev^2 / sum(pc1$sdev^2)

ggplot(pcr, aes(PC1, PC2, color=group)) + geom_point(size=10, alpha=0.6, shape = 19) + theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12),
        legend.position = "right") + 
  labs(
    x = paste0("PC1 (", round(pc_var[1] * 100, 2), "% variance)"),
    y = paste0("PC2 (", round(pc_var[2] * 100, 2), "% variance)"),
    title = "PCA Plot"
  )

ggplot(pcr1, aes(PC1, PC2, color=group)) + geom_point(size=10, alpha=0.6, shape = 19) + theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12),
        legend.position = "right") + 
  labs(
    x = paste0("PC1 (", round(pc1_var[1] * 100, 2), "% variance)"),
    y = paste0("PC2 (", round(pc1_var[2] * 100, 2), "% variance)"),
    title = "PCA Plot"
  )


########################### used with Data 2 that are based on Anova analysis ###############

data2 <- read_excel("/Users/mortezaabyadeh/Desktop/Transcript.xlsx")
data2 <- as.data.frame(data2)

head(data2)
dim(data2)
summary(data2)
data <- subset(data2, ANOVA < 0.05)


data[,1] <- rownames(data)
data <- data[,-1]
dim(data)
data <- data[,-13]
head(data)
data <- na.omit(data)
data <- data.frame(lapply(data, function(x) if(is.numeric(x)) as.integer(x) else x))
head(data)
summary(data)
dim(data)
gr <- c(rep("CTRL", 3), rep("Eto", 3), rep("FTY", 3), rep("FTY-Eto", 3))
gr <- as.factor(gr)

dim(data)
head(data, 5)
data <- na.omit(data)
head(data)

colnames(data) <- gr
colData <- data.frame(Group = gr)
head(data)

rownames(colData) <- colnames(data)

cds <- DESeqDataSetFromMatrix(countData = data, colData = colData, design = ~ Group)
cds
cds <- DESeq(cds)


data.norm <- log2(1+counts(cds, normalized=T))
head(data.norm)
boxplot(data.norm)

data.mean.center <- t(scale(t(data.norm), scale = F))

### PCA

pc <- prcomp(data.mean.center)
pc1 <- prcomp(data.norm)
plot(pc1)
plot(pc)
head(pc$rotation)
pcr <- data.frame(pc$r)
pcr1 <- data.frame(pc1$r)
head(pcr)
head(data)
pcr1$group <- gr
pcr$group <- gr
head(pcr)


pc_var <- pc$sdev^2 / sum(pc$sdev^2)
pc_var

pc1_var <- pc1$sdev^2 / sum(pc1$sdev^2)

ggplot(pcr, aes(PC1, PC2, color=group)) + geom_point(size=10, alpha=0.6, shape = 19) + theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12),
        legend.position = "right") + 
  labs(
    x = paste0("PC1 (", round(pc_var[1] * 100, 2), "% variance)"),
    y = paste0("PC2 (", round(pc_var[2] * 100, 2), "% variance)"),
    title = "PCA Plot"
  )

ggplot(pcr1, aes(PC1, PC2, color=group)) + geom_point(size=10, alpha=0.6, shape = 19) + theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12),
        legend.position = "right") + 
  labs(
    x = paste0("PC1 (", round(pc1_var[1] * 100, 2), "% variance)"),
    y = paste0("PC2 (", round(pc1_var[2] * 100, 2), "% variance)"),
    title = "PCA Plot"
  )


####### It seems it was better before anova! 




############## heatmap using the anova p-value ###################


data2 <- read_excel("/Users/mortezaabyadeh/Desktop/Transcript.xlsx")
data2 <- as.data.frame(data2)

head(data2)
dim(data2)



colnames(data2)[11:13] <- sub("-", "_", colnames(data2)[11:13])

head(data2)









heatmap_maker <- function(df){
  df <- df
  df_filtered <- df[df$`ANOVA` < 0.05, ]
  df_filtered <- df_filtered[complete.cases(df_filtered), ]
  i_subset <- df_filtered[,2:13]
  i_matrix <- as.matrix(i_subset)
  i_log <- log2(i_matrix + 1)
  i <- t(scale(t(i_log), scale = FALSE))
  i_annotation_col <- data.frame(Group= c(sub("-.*", "", colnames(i))))
  rownames(i_annotation_col) <- colnames(i)
  i_head <- pheatmap(i,
                     cluster_cols = TRUE,  
                     scale = "row",  
                     color = colorRampPalette(c("green", "black", "red"))(100), 
                     main = "Heatmap of Normalized Counts-Transcriptome",  
                     fontsize = 15,  
                     fontsize_row = 12,  
                     fontsize_col = 17,  
                     border_color = NA,  
                     angle_col = 45,
                     annotation_col = i_annotation_col,
                     show_rownames = FALSE,
                     width = 7,
                     height = 7
  )
}

dim(data)
#data <- data[,-13]
head(data)
heatmap_maker(data2)

correlation_maker <- function(df) {
  df1 <- df
  df_filtered1 <- df1[df1$`ANOVA` < 0.05, ]
  df_filtered1 <- df_filtered1[complete.cases(df_filtered1), ]
  i_subset1 <- df_filtered1[,2:13]
  i_matrix1 <- as.matrix(i_subset1)
  i_log1 <- log2(i_matrix1 + 1)
  i1 <- t(scale(t(i_log1), scale = FALSE))
  
  corr_matrix <- cor(i1)
  corr_melted <- melt(corr_matrix)
  
  ggplot(data = corr_melted, aes(Var1, Var2, fill = value)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0, limit = c(-1, 1), space = "Lab", 
                         name = "Correlation") +
    theme_minimal() + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1),
          axis.text.y = element_text(size = 12)) +
    coord_fixed() +
    ggtitle("Correlation") +
    labs(x = "", y = "")
}
correlation_maker(data2)



head(data2)


correlation_makers <- function(df) {
  df1 <- df
  df_filtered1 <- df1[df1$ANOVA < 0.05, ]
  df_filtered1 <- df_filtered1[complete.cases(df_filtered1), ]
  i_subset1 <- df_filtered1[, 2:13]  
  i_matrix1 <- as.matrix(i_subset1)
  i_log1 <- log2(i_matrix1 + 1)
  i1 <- t(scale(t(i_log1), scale = FALSE))
  
  corr_matrix <- cor(i1)
  
  corr_matrix[upper.tri(corr_matrix, diag = FALSE)] <- NA 
  
  corr_melted <- melt(corr_matrix, na.rm = TRUE)
  
  ggplot(data = corr_melted, aes(Var1, Var2, fill = value)) +
    geom_tile(color = "white", na.rm = TRUE) +
    geom_text(aes(label = round(value, 2)), size = 4, na.rm = TRUE) +  
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0, limit = c(-1, 1), space = "Lab", 
                         name = "Correlation") +
    theme_minimal() + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 14, hjust = 1),
          axis.text.y = element_text(size = 14),
          panel.grid = element_blank()) +  
    scale_y_discrete(position = "right") +  
    coord_fixed() +
    ggtitle("Transcriptome Correlation Heatmap") +
    labs(x = "", y = "")
}




correlation_makers(data2)



############################ Proteome data ########################

proteome <- read_excel("/Users/mortezaabyadeh/Desktop/proteome_Anova_Significant.xlsx")

head(proteome)
data <- read_excel("proteome_Anova_Significant.xlsx")

class(proteome)
proteome <- as.data.frame(proteome)


p <- proteome[, 1:13]
colnames(p)


colnames(p) <- c("Gene","CTRL-1","CTRL-2", "CTRL-3","Eto-1","Eto-2","Eto-3","FTY-1","FTY-2","FTY-3","FTY_Eto-1","FTY_Eto-2","FTY_Eto-3")

colnames(p)


heatmap_maker1 <- function(df){
  df_filtered <- df
  df_filtered <- df_filtered[complete.cases(df_filtered), ]
  i_subset <- df_filtered[,2:13]
  i_matrix <- as.matrix(i_subset)
  i_log <- log2(i_matrix + 1)
  i <- t(scale(t(i_log), scale = FALSE))
  i_annotation_col <- data.frame(Group= c(sub("-.*", "", colnames(i))))
  rownames(i_annotation_col) <- colnames(i)
  i_head <- pheatmap(i,
                     cluster_cols = FALSE,  
                     scale = "row",  
                     color = colorRampPalette(c("green", "black", "red"))(100), 
                     main = "Heatmap of Normalized Counts-proteome",  
                     fontsize = 15,  
                     fontsize_row = 12,  
                     fontsize_col = 17,  
                     border_color = NA,  
                     angle_col = 45,
                     annotation_col = i_annotation_col,
                     show_rownames = FALSE,
                     width = 7,
                     height = 7
  )
}
heatmap_maker1(p)



correlation_makers1 <- function(df) {
  df_filtered1 <- df
  df_filtered1 <- df_filtered1[complete.cases(df_filtered1), ]
  i_subset1 <- df_filtered1[, 2:13]  
  i_matrix1 <- as.matrix(i_subset1)
  i_log1 <- log2(i_matrix1 + 1)
  i1 <- t(scale(t(i_log1), scale = FALSE))
  
  corr_matrix <- cor(i1)
  
  corr_matrix[upper.tri(corr_matrix, diag = FALSE)] <- NA 
  
  corr_melted <- melt(corr_matrix, na.rm = TRUE)
  
  ggplot(data = corr_melted, aes(Var1, Var2, fill = value)) +
    geom_tile(color = "white", na.rm = TRUE) +
    geom_text(aes(label = round(value, 2)), size = 4, na.rm = TRUE) +  
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0, limit = c(-1, 1), space = "Lab", 
                         name = "Correlation") +
    theme_minimal() + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 14, hjust = 1),
          axis.text.y = element_text(size = 14),
          panel.grid = element_blank()) +  
    scale_y_discrete(position = "right") +  
    coord_fixed() +
    ggtitle("Proteome Correlation Heatmap") +
    labs(x = "", y = "")
}

correlation_makers1(p)


