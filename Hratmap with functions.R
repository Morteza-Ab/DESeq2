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

#data_proteome <- read_excel("set5.xlsx")
data_proteome <- read_excel("set5.xlsx", sheet = "DAPs")



head(data_proteome)
sorted_pro <- data_proteome %>% arrange(data_proteome$log2FC)
head(sorted_pro)
down_pr <- sorted_pro[ 1:10, ]
head(down_pr,13)
up_pr <- tail(sorted_pro,10)
head(up_pr,13)
data_pr <- rbind(down_pr, up_pr)
head(data_pr,25)
data_pr <- data_pr[,1:7]
head(data_pr)

data_pr <- as.data.frame(data_pr)
rownames(data_pr) <- data_pr$Gene

data_pr[,2:7] <- lapply(data_pr[,2:7], as.numeric)
data_pr <- data_pr[,-1]
head(data_pr)

pr_matrix <- as.matrix(data_pr)

print(dim(pr_matrix))  
print(colnames(pr_matrix))  

pr_matrix <- log2(pr_matrix + 1)

pr_data.mean.center <- t(scale(t(pr_matrix), scale = F))

annotation_col <- data.frame(Group = c("Eto", "Eto", "Eto", "FTY_Eto", "FTY_Eto", "FTY_Eto"))
rownames(annotation_col) <- colnames(pr_matrix)

print(dim(annotation_col))

pheatmap(pr_matrix,
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

pheatmap(pr_data.mean.center,
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


### for genecard-pr
df_pr <- read.csv("p_genecard.csv",row.names="Gene")

df_pr <- as.data.frame(df_pr)
head(df_pr,20)

dim(df_pr)

pheatmap(pr_data.mean.center,
         cluster_cols = TRUE,  
         cluster_rows = TRUE,  
         scale = "row",  
         color = colorRampPalette(c("blue", "white", "red"))(100), 
         main = "Heatmap of Normalized Counts for top 20 DAPs",  
         fontsize = 13,  
         fontsize_row = 8,  # Adjust fontsize_row to ensure visibility
         fontsize_col = 13,  
         display_numbers = TRUE,  
         number_color = "black",
         border_color = NA,  
         angle_col = 45,
         rownames_force = TRUE,  # Force display of row names
         annotation_col = annotation_col,
         annotation_row = df_pr
)

##save by 8 * 20



### for uniprot-pr

df_pr1 <- read.csv("p_uni.csv",row.names="Gene")

df_pr1 <- as.data.frame(df_pr1)
head(df_pr1,20)

dim(df_pr1)

pheatmap(pr_data.mean.center,
         cluster_cols = TRUE,  
         cluster_rows = TRUE,  
         scale = "row",  
         color = colorRampPalette(c("blue", "white", "red"))(100), 
         main = "Heatmap of Normalized Counts for top 20 DAPs",  
         fontsize = 13,  
         fontsize_row = 13,  # Adjust fontsize_row to ensure visibility
         fontsize_col = 13,  
         display_numbers = TRUE,  
         number_color = "black",
         border_color = NA,  
         angle_col = 45,
         rownames_force = TRUE,  # Force display of row names
         annotation_col = annotation_col,
         annotation_row = df_pr1
)

# save 10in * 23






########### Transcriptome data

data_transc <- read_excel("clean_transcriptome_FTY_Eto_vs_Eto.xlsx", sheet = "DEGs")



head(data_transc)
sorted_tran <- data_transc %>% arrange(data_transc$log2FC)
head(sorted_tran)
down_tr <- sorted_tran[ 1:10, ]
head(down_tr,13)
up_tr <- tail(sorted_tran,10)
head(up_tr,13)
data_tr <- rbind(down_tr, up_tr)
head(data_tr,25)
data_tr <- data_tr[,1:7]
head(data_tr)

data_tr <- as.data.frame(data_tr)
rownames(data_tr) <- data_tr$Gene

data_tr[,2:7] <- lapply(data_tr[,2:7], as.numeric)
data_tr <- data_tr[,-1]
head(data_tr)

tr_matrix <- as.matrix(data_tr)

print(dim(tr_matrix))  
print(colnames(tr_matrix))  

tr_matrix <- log2(tr_matrix + 1)

tr_data.mean.center <- t(scale(t(tr_matrix), scale = F))

annotation_col <- data.frame(Group = c("Eto", "Eto", "Eto", "FTY_Eto", "FTY_Eto", "FTY_Eto"))
rownames(annotation_col) <- colnames(tr_matrix)

print(dim(annotation_col))




### for genecard-tr
df_tr <- read.csv("t_genecard.csv",row.names="Gene")

df_tr <- as.data.frame(df_tr)
head(df_tr,20)

dim(df_tr)

pheatmap(tr_data.mean.center,
         cluster_cols = TRUE,  
         cluster_rows = TRUE,  
         scale = "row",  
         color = colorRampPalette(c("blue", "white", "red"))(100), 
         main = "Heatmap of Normalized Counts for top 20 DEGs",  
         fontsize = 13,  
         fontsize_row = 8,  # Adjust fontsize_row to ensure visibility
         fontsize_col = 13,  
         display_numbers = TRUE,  
         number_color = "black",
         border_color = NA,  
         angle_col = 45,
         rownames_force = TRUE,  # Force display of row names
         annotation_col = annotation_col,
         annotation_row = df_tr
)

##save by 8 * 20



### for uniprot-tr

df_tr1 <- read.csv("t_uni.csv",row.names="Gene")

df_tr1 <- as.data.frame(df_tr1)
head(df_tr1,20)

dim(df_tr1)

pheatmap(tr_data.mean.center,
         cluster_cols = TRUE,  
         cluster_rows = TRUE,  
         scale = "row",  
         color = colorRampPalette(c("blue", "white", "red"))(100), 
         main = "Heatmap of Normalized Counts for top 20 DEGs",  
         fontsize = 13,  
         fontsize_row = 13,  # Adjust fontsize_row to ensure visibility
         fontsize_col = 13,  
         display_numbers = TRUE,  
         number_color = "black",
         border_color = NA,  
         angle_col = 45,
         rownames_force = TRUE,  # Force display of row names
         annotation_col = annotation_col,
         annotation_row = df_tr1
)

# save by 8 *20







