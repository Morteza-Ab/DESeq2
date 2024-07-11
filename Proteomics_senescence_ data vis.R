setwd("/Users/mortezaabyadeh/Desktop")
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


data <- read_excel("proteome_Anova_Significant.xlsx")

class(data)
colnames(data)

data <- data[, 2:14]
colnames(data)
data[,1] <- rownames(data)
data <- data[,-1]

data <- na.omit(data)
data <- data.frame(lapply(data, function(x) if(is.numeric(x)) as.integer(x) else x))

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
plot(pc)
head(pc$rotation)
pcr <- data.frame(pc$r)
head(pcr)
head(data)
pcr$group <- gr
head(pcr)


pc_var <- pc$sdev^2 / sum(pc$sdev^2)
pc_var

ggplot(pcr, aes(PC1, PC2, color=group)) + geom_point(size=10, alpha=0.6, shape = 19) + theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12),
        legend.position = "right") + 
  labs(
    x = paste0("PC1 (", round(pca_var[1] * 100, 2), "% variance)"),
    y = paste0("PC2 (", round(pca_var[2] * 100, 2), "% variance)"),
    title = "PCA Plot"
  )


### 3-D PCA
fig <- plot_ly(data = pcr, x = ~PC1, y = ~PC2, z = ~PC3, type = "scatter3d", mode = "markers", color = ~group) %>%
  layout(
    scene = list(
      xaxis = list(title = paste0("PC1 (", round(pc_var[1] * 100, 2), "% variance)")),
      yaxis = list(title = paste0("PC2 (", round(pc_var[2] * 100, 2), "% variance)")),
      zaxis = list(title = paste0("PC3 (", round(pc_var[3] * 100, 2), "% variance)"))
    ),
    title = "3D PCA Plot"
  )


# I have saved the 3-D pca and send the html file just need to deploy it soon!
saveWidget(fig, file = "/Users/mortezaabyadeh/Desktop/3d_PCA_proteomics.html")

### Correlation
pheatmap(data.mean.center)
?pheatmap
pheatmap(data.norm)
pheatmap(cor(data.mean.center))
pheatmap(cor(data.mean.center))
head(data.mean.center)


### Heatmap only for the DAPs that already obtained from ANOVA and related posthoc!


excel_file <- "/Users/mortezaabyadeh/Desktop/merged_posthoc_proteome_anova.xlsx"

datasets <- lapply(excel_sheets(excel_file), function(sheet) read_excel(excel_file, sheet = sheet))

names(datasets) <- excel_sheets(excel_file)


for (sheet_name in names(datasets)) {
  cat("Dataset from sheet:", sheet_name, "\n")
  print(head(datasets[[sheet_name]]))  
}

colnames(datasets$'Eto-ctrl')

control_cols = c("Abundances (Normalized): F1: Sample, Control", "Abundances (Normalized): F2: Sample, Control", "Abundances (Normalized): F3: Sample, Control")
etoposide_cols = c("Abundances (Normalized): F4: Sample, Etoposide", "Abundances (Normalized): F5: Sample, Etoposide", "Abundances (Normalized): F6: Sample, Etoposide")
fty_cols = c("Abundances (Normalized): F10: Sample, FTY720", "Abundances (Normalized): F11: Sample, FTY720", "Abundances (Normalized): F12: Sample, FTY720")
fty_eto_cols = c("Abundances (Normalized): F7: Sample, FTY720_plus_Etoposide", "Abundances (Normalized): F8: Sample, FTY720_plus_Etoposide","Abundances (Normalized): F9: Sample, FTY720_plus_Etoposide")

datasets$'Eto-ctrl' <- datasets$'Eto-ctrl'[ , c("Gene", "p-adj", "FC_Eto/CTRL_x", control_cols, etoposide_cols, "ANOVA")]

colnames(datasets$'Eto-ctrl')
new_colnames <- c("Gene", "Adj_P_Value", "Fold_Change_Eto_CTRL", "CTRL-1", "CTRL-2", "CTRL-3", "Eto-1", "Eto-2", "Eto-3", "ANOVA")
colnames(datasets$'Eto-ctrl') <- new_colnames
data3 <- datasets$'Eto-ctrl'[,c("Gene", "CTRL-1", "CTRL-2", "CTRL-3", "Eto-1", "Eto-2", "Eto-3")]
dim(data3)


# I removed the duplicated within gene names
data3_unique <- data3 %>%
  distinct(Gene, .keep_all = TRUE)
dim(data3_unique)

data3_unique[c("CTRL-1", "CTRL-2", "CTRL-3")] <- lapply(data3_unique[c("CTRL-1", "CTRL-2", "CTRL-3")], as.numeric)
data3_unique[c("Eto-1", "Eto-2", "Eto-3")] <- lapply(data3_unique[c("Eto-1", "Eto-2", "Eto-3")], as.numeric)

rownames(data3_unique) <- data3_unique$Gene

data3_unique <- data3_unique[,c("CTRL-1", "CTRL-2", "CTRL-3", "Eto-1", "Eto-2", "Eto-3")]

data_matrix <- as.matrix(data3_unique)

print(dim(data_matrix))  
print(colnames(data_matrix)) 


data_matrix <- log2(data_matrix + 1)

data.mean.center <- t(scale(t(data_matrix), scale = F))



annotation_col <- data.frame(Group = c("CTRL", "CTRL", "CTRL", "Eto", "Eto", "Eto"))
rownames(annotation_col) <- colnames(data_matrix)
print(dim(annotation_col))

pheatmap(data_matrix,
         cluster_cols = TRUE,  
         cluster_rows = TRUE,  
         scale = "row",  
         color = colorRampPalette(c("blue", "white", "red"))(100), 
         main = "Heatmap of Normalized Counts (Eto/CTRL)",  
         fontsize = 10,  
         fontsize_row = 8,  
         fontsize_col = 8,  
         #display_numbers = TRUE,  
         #number_color = "black",
         border_color = NA,  
         angle_col = 45,
         #rownames_force = TRUE, 
         annotation_col = annotation_col
)


##################### make only for first 10 genes
sorted_data <- datasets$'Eto-ctrl' %>%
  arrange("Adj_P_Value")
colnames(datasets$'Eto-ctrl')
top_data <- head(sorted_data, 10)

data4 <- top_data[,c("Gene","CTRL-1", "CTRL-2", "CTRL-3", "Eto-1", "Eto-2", "Eto-3")]
dim(data4)
data4_unique <- data4 %>%
  distinct(Gene, .keep_all = TRUE)
dim(data4_unique)

data4_unique <- as.data.frame(data4_unique)
rownames(data4_unique) <- data4_unique$Gene
head(data4_unique)
data4_unique[c("CTRL-1", "CTRL-2", "CTRL-3")] <- lapply(data4_unique[c("CTRL-1", "CTRL-2", "CTRL-3")], as.numeric)
data4_unique[c("Eto-1", "Eto-2", "Eto-3")] <- lapply(data4_unique[c("Eto-1", "Eto-2", "Eto-3")], as.numeric)

data4_unique <- data4_unique[,-1]
head(data4_unique,10)

data_matrix <- as.matrix(data4_unique)

print(dim(data_matrix))  
print(colnames(data_matrix))  

data_matrix <- log2(data_matrix + 1)

data.mean.center <- t(scale(t(data_matrix), scale = F))

annotation_col <- data.frame(Group = c("CTRL", "CTRL", "CTRL", "Eto", "Eto", "Eto"))
rownames(annotation_col) <- colnames(data_matrix)
print(dim(annotation_col))

pheatmap(data_matrix,
         cluster_cols = TRUE,  
         cluster_rows = TRUE,  
         scale = "row",  
         color = colorRampPalette(c("blue", "white", "red"))(100), 
         main = "Heatmap of Normalized Counts for top 10 DAPs",  
         fontsize = 13,  
         fontsize_row = 11,  
         fontsize_col = 13,  
         display_numbers = TRUE,  
         number_color = "black",
         border_color = NA,  
         angle_col = 45,
         rownames_force = TRUE,  # Force display of row names
         annotation_col = annotation_col
)


### adding function of proteins to heatmap
df <- read.csv("Function.csv", row.names="Gene")
head(df)

df <- as.data.frame(df)

pheatmap(data_matrix,
         cluster_cols = TRUE,  
         cluster_rows = TRUE,  
         scale = "row",  
         color = colorRampPalette(c("blue", "white", "red"))(100), 
         main = "Heatmap of Normalized Counts for top 10 DAPs",  
         fontsize = 13,  
         fontsize_row = 11,  
         fontsize_col = 13,  
         display_numbers = TRUE,  
         number_color = "black",
         border_color = NA,  
         angle_col = 45,
         rownames_force = TRUE,  # Force display of row names
         annotation_col = annotation_col,
         annotation_row = df
)



### volcano

df <- "/Users/mortezaabyadeh/Desktop/posthoc_results2.xlsx"

datasets1 <- lapply(excel_sheets(df), function(sheet) read_excel(df, sheet = sheet))

names(datasets1) <- excel_sheets(df)  # Assign sheet names to list names


for (sheet_name in names(datasets1)) {
  cat("Dataset from sheet:", sheet_name, "\n")
  print(head(datasets1[[sheet_name]]))  
}

df1 <- datasets1$"Control_vs_Etoposide"[,c("Gene", "p-adj", "log2_FC_Eto/CTRL")]


df1$significance <- ifelse(df1$`p-adj` < 0.05 & df1$`log2_FC_Eto/CTRL` > 0, "Significant Up",
                           ifelse(df1$`p-adj` < 0.05 & df1$`log2_FC_Eto/CTRL` < 0, "Significant Down", "Not Significant"))

ggplot(df1, aes(x = `log2_FC_Eto/CTRL`, y = -log10(`p-adj`))) +
  geom_point(aes(
    color = significance,
    shape = significance
  ), size = 2, alpha = 0.6, show.legend = c(size = FALSE)) +
  scale_color_manual(values = c("Significant Up" = "red", "Significant Down" = "blue", "Not Significant" = "grey")) +
  scale_shape_manual(values = c("Significant Up" = 20, "Significant Down" = 20, "Not Significant" = 19)) +
  labs(
    x = "log2(Fold Change)",
    y = "-log10(P-value)",
    title = "Volcano Plot"
  ) +
  scale_y_continuous(limits = c(0, 4), expand = c(0, 0)) +
  theme_minimal() +
  theme(
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 18)
  )





#### Volcano
df1$diffexpressed <- ifelse(df1$`p-adj` < 0.05 & df1$`log2_FC_Eto/CTRL` > 0, "Up",
                           ifelse(df1$`p-adj` < 0.05 & df1$`log2_FC_Eto/CTRL` < 0, "Down", "No"))

df1$delabel <- df1$Gene
head(df1)
write.csv(df1, file = "df1_data.csv", row.names = FALSE)

df <- read.csv("df1_data.csv")

head(df)



ggplot(data = df, aes(x = log2_FC_Eto.CTRL, y = -log10(p.adj), col = diffexpressed)) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
  geom_point(size = 2) + 
  scale_color_manual(values = c("#00AFBB", "grey", "#bb0c00"), # to set the colours of our variable  
                     labels = c("Downregulated", "Not significant", "Upregulated")) + # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)
  coord_cartesian(ylim = c(0, 5), xlim = c(-10, 10)) + # since some genes can have minuslog10padj of inf, we set these limits
  #labs(color = 'Severe', #legend_title, 
       #x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) + 
  scale_x_continuous(breaks = seq(-10, 10, 2)) + # to customise the breaks in the x axis
  ggtitle('Proteome Eto vs. CTRL') # Plot title 
  #geom_text_repel(max.overlaps = Inf) # To show all labels 

