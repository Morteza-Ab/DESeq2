library(DESeq2)
library(readxl)
library(dplyr)
library(openxlsx)
library(pheatmap)
library(data.table)
library(ggplot2)
library(DESeq2)
library(openxlsx) # not needed here but kept for other analysis
library(readxl)
library(plotly)
library(pheatmap)
library(htmlwidgets)
library(dplyr)
library(RColorBrewer) 
library(ggrepel)
library(tidyverse)
library(reshape2)
library(shiny)
library(pheatmap)
install.packages("ggalt")
library(ggalt)

setwd("/Users/mortezaabyadeh/Desktop")

df <- read_excel("clean_transcriptome_Eto_vs_CTRL.xlsx")
head(df,2)

df_filtered <- df[df$p_adj < 0.05, ]
numeric_columns <- df_filtered[,2:7]
numeric_columns <- data.frame(lapply(numeric_columns, function(x) if(is.numeric(x)) as.integer(x) else x))
numeric_columns <- na.omit(numeric_columns)
gr <- sub("\\-.*", "", colnames(numeric_columns))
gr <- as.factor(gr)
colnames(numeric_columns) <- gr
col_data <- data.frame(Group = gr)
cds <- DESeqDataSetFromMatrix(countData = numeric_columns, colData = col_data, design = ~ Group)
data.norm <- log2(1 + counts(cds, normalized = F))
pc <- prcomp(t(data.norm))
pcr <- data.frame(pc$x)
pcr$group <- gr
pc_var <- pc$sdev^2 / sum(pc$sdev^2)
pcr

pca1 <- read_excel("pca.xls")
ggplot(pca1, aes(PC1, PC2, color = group)) +
  geom_point(size = 9, alpha = 0.6) +
  theme_bw() +
  labs(x = paste0("PC1 (", round(pc_var[1] * 100, 2), "% variance)"),
       y = paste0("PC2 (", round(pc_var[2] * 100, 2), "% variance)"),
       title = "PCA Plot")





plot_ly(data = pca1, x = ~PC1, y = ~PC2, z = ~PC3, type = "scatter3d", mode = "markers", color = ~group) %>%
    layout(
      scene = list(
        xaxis = list(title = paste0("PC1 (", round(pc_var[1] * 100, 2), "% variance)")),
        yaxis = list(title = paste0("PC2 (", round(pc_var[2] * 100, 2), "% variance)")),
        zaxis = list(title = paste0("PC3 (", round(pc_var[3] * 100, 2), "% variance)"))
      ),
      title = "3D PCA Plot"
)


ggplot(pca1, aes(PC1, PC2, color = group)) +
  geom_point(size = 9, alpha = 0.6) +
  geom_encircle(aes(group = group), color = "black", size = 1, expand = 0.1) +
  theme_bw() +
  labs(x = paste0("PC1 (", round(pc_var[1] * 100, 2), "% variance)"),
       y = paste0("PC2 (", round(pc_var[2] * 100, 2), "% variance)"),
       title = "PCA Plot")




light_colors <- c("#FFB3B3", "#B3D9FF", "#B3FFB3", "#FFD9B3", "#FFB3FF", "#FFFFB3")

ggplot(pca1, aes(PC1, PC2, color = group)) +
  geom_point(size = 9, alpha = 0.6) +
  geom_encircle(aes(PC1, PC2, group = group, fill = group), 
                color = NA,  # Remove the color argument here
                size = 1.5, 
                expand = 0.1,
                alpha = 0.4) +  # Adjust transparency for a lighter effect
  theme_bw() +
  labs(x = paste0("PC1 (", round(pc_var[1] * 100, 2), "% variance)"),
       y = paste0("PC2 (", round(pc_var[2] * 100, 2), "% variance)"),
       title = "PCA Plot") +
  scale_fill_manual(values = light_colors) +
  scale_color_manual(values = light_colors)
