library(shiny)
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