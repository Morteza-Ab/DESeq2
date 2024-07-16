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
library(reshape2)


file_path <- "/Users/mortezaabyadeh/Desktop/merged_posthoc_proteome_anova.xlsx"
sheet_names <- excel_sheets(file_path)
datasets <- lapply(sheet_names, function(sheet) read_excel(file_path, sheet = sheet))
names(datasets) <- sheet_names

lapply(datasets, head)
colnames(datasets$`Eto-ctrl`)

datasets$'Eto-ctrl'["log2_FC_Eto/CTRL_x"] <- log2(datasets$'Eto-ctrl'["FC_Eto/CTRL_x"])
head(datasets$'Eto-ctrl'[, c('FC_Eto/CTRL_x', 'log2_FC_Eto/CTRL_x')], 2)

datasets$'FTY-CTRL'["log2_FC_FTY/CTRL_x"] <- log2(datasets$"FTY-CTRL"["FC_FTY/CTRL_x"])
datasets$'FTY-Eto-CTRL'["log2_FC_FTY_Eto/CTRL_x"] <- log2(datasets$'FTY-Eto-CTRL'["FC_FTY_Eto/CTRL_x"])
datasets$'FTY-Eto'["log2_FC_FTY/Eto_x"] <- log2(datasets$'FTY-Eto'["FC_FTY/Eto_x"])
datasets$'FTY-Eto-Eto'["log2_FC_FTY_Eto/Eto_x"] <- log2(datasets$'FTY-Eto-Eto'["FC_FTY_Eto/Eto_x"])
datasets$'FTY-Eto-FTY'["log2_FC_FTY_Eto/FTY_x"] <- log2(datasets$'FTY-Eto-FTY'["FC_FTY_Eto/FTY_x"])


control_cols = c("Abundances (Normalized): F1: Sample, Control", "Abundances (Normalized): F2: Sample, Control", "Abundances (Normalized): F3: Sample, Control")
etoposide_cols = c("Abundances (Normalized): F4: Sample, Etoposide", "Abundances (Normalized): F5: Sample, Etoposide", "Abundances (Normalized): F6: Sample, Etoposide")
fty_cols = c("Abundances (Normalized): F10: Sample, FTY720", "Abundances (Normalized): F11: Sample, FTY720", "Abundances (Normalized): F12: Sample, FTY720")
fty_eto_cols = c("Abundances (Normalized): F7: Sample, FTY720_plus_Etoposide", "Abundances (Normalized): F8: Sample, FTY720_plus_Etoposide","Abundances (Normalized): F9: Sample, FTY720_plus_Etoposide")

datasets$'Eto-ctrl' <- datasets$'Eto-ctrl'[ , c("Gene", "p-adj", "FC_Eto/CTRL_x", "log2_FC_Eto/CTRL_x" , control_cols, etoposide_cols, "ANOVA")]
datasets$"FTY-CTRL" <- datasets$'FTY-CTRL'[ , c("Gene", "p-adj", "FC_FTY/CTRL_x", "log2_FC_FTY/CTRL_x",control_cols, fty_cols, "ANOVA")]
datasets$'FTY-Eto-CTRL' <- datasets$'FTY-Eto-CTRL'[,c("Gene", "p-adj", "FC_FTY_Eto/CTRL_x", "log2_FC_FTY_Eto/CTRL_x",control_cols, fty_eto_cols, "ANOVA")]
datasets$"FTY-Eto" <- datasets$"FTY-Eto"[,c("Gene", "p-adj", "FC_FTY/Eto_x", "log2_FC_FTY/Eto_x",fty_cols, etoposide_cols, "ANOVA")]
datasets$"FTY-Eto-Eto" <- datasets$"FTY-Eto-Eto"[,c("Gene", "p-adj", "FC_FTY_Eto/Eto_x", "log2_FC_FTY_Eto/Eto_x",fty_eto_cols, etoposide_cols, "ANOVA")]
datasets$'FTY-Eto-FTY' <- datasets$'FTY-Eto-FTY'[,c("Gene", "p-adj", "FC_FTY_Eto/FTY_x", "log2_FC_FTY_Eto/FTY_x",fty_eto_cols, fty_cols, "ANOVA")]


new_colnames_Eto_ctrl <- c("Gene", "Adj_P_Value", "Fold_Change_Eto_CTRL", "log2Fold_change_Eto_CTRL", "CTRL-1", "CTRL-2", "CTRL-3", "Eto-1", "Eto-2", "Eto-3", "ANOVA")
new_colnames_FTY_CTRL <- c("Gene", "Adj_P_Value", "Fold_Change_FTY_CTRL", "log2Fold_change_FTY_CTRL", "CTRL-1", "CTRL-2", "CTRL-3", "FTY-1", "FTY-2", "FTY-3", "ANOVA")
new_colnames_FTY_Eto_CTRL <- c("Gene", "Adj_P_Value", "Fold_Change_FTY_Eto_CTRL", "log2Fold_change_FTY_Eto_CTRL", "CTRL-1", "CTRL-2", "CTRL-3", "FTY_Eto-1", "FTY_Eto-2", "FTY_Eto-3", "ANOVA")
new_colnames_FTY_Eto <- c("Gene", "Adj_P_Value", "Fold_Change_FTY_Eto", "log2Fold_change_FTY_Eto", "FTY-1", "FTY-2", "FTY-3", "Eto-1", "Eto-2", "Eto-3", "ANOVA")
new_colnames_FTY_Eto_Eto <- c("Gene", "Adj_P_Value", "Fold_Change_FTY_Eto_Eto", "log2Fold_change_FTY_Eto_Eto", "FTY_Eto-1", "FTY_Eto-2", "FTY_Eto-3", "Eto-1", "Eto-2", "Eto-3", "ANOVA")
new_colnames_FTY_Eto_FTY <- c("Gene", "Adj_P_Value", "Fold_Change_FTY_Eto_FTY", "log2Fold_change_FTY_Eto_FTY", "FTY_Eto-1", "FTY_Eto-2", "FTY_Eto-3", "FTY-1", "FTY-2", "FTY-3", "ANOVA")

colnames(datasets$'Eto-ctrl') <- new_colnames_Eto_ctrl
colnames(datasets$'FTY-CTRL') <- new_colnames_FTY_CTRL
colnames(datasets$'FTY-Eto-CTRL') <- new_colnames_FTY_Eto_CTRL
colnames(datasets$'FTY-Eto') <- new_colnames_FTY_Eto
colnames(datasets$'FTY-Eto-Eto') <- new_colnames_FTY_Eto_Eto
colnames(datasets$'FTY-Eto-FTY') <- new_colnames_FTY_Eto_FTY

colnames(datasets$'Eto-ctrl')
colnames(datasets$'FTY-CTRL')
colnames(datasets$'FTY-Eto-CTRL')
colnames(datasets$'FTY-Eto')
colnames(datasets$'FTY-Eto-Eto')
colnames(datasets$'FTY-Eto-FTY')


dataset1 <- datasets$'Eto-ctrl' %>%
  arrange(log2Fold_change_Eto_CTRL)
head(datasets$`Eto-ctrl`,2)
head(dataset1)
tail(dataset1)
colnames(dataset1)
dataset1 <- datasets$'Eto-ctrl' %>% arrange(log2Fold_change_Eto_CTRL)
dataset2 <- datasets$`FTY-CTRL`%>% arrange(log2Fold_change_FTY_CTRL)
dataset3 <- datasets$'FTY-Eto-CTRL' %>% arrange(log2Fold_change_FTY_Eto_CTRL)
dataset4 <- datasets$'FTY-Eto' %>% arrange(log2Fold_change_FTY_Eto)
dataset5 <- datasets$'FTY-Eto-Eto' %>% arrange(log2Fold_change_FTY_Eto_Eto)
dataset6 <- datasets$'FTY-Eto-FTY' %>% arrange(log2Fold_change_FTY_Eto_FTY)

datasets <- list(dataset1, dataset2, dataset3, dataset4, dataset5, dataset6)


sapply(datasets, function(df) class(df$Gene))
gsn_present <- sapply(datasets, function(df) any(df$Gene == "GSN"))
print(gsn_present)
sum(sapply(datasets, function(df) sum(df$Gene == "GSN")))
datasets <- lapply(datasets, function(df) df[df$Gene != "GSN", ])
sum(sapply(datasets, function(df) sum(df$Gene == "GSN")))
print_dimensions(datasets)
head(dataset1)

for (i in seq_along(datasets)) {
  write.csv(datasets[[i]], file = paste0("cleaned", i, ".csv"), row.names = FALSE)
}
######################################################### Data has been saved /
# will be loaded a few lines below
dataset1 <- read.delim("cleaned1.csv", header = TRUE)
head(dataset1)
rownames(dataset1) <- dataset1[,1]
head(dataset1[,1])
dataset1 <- dataset1[,-1]
rownames(dataset1)

print_dimensions <- function(dim_list) {
  for (i in dim_list) {
    print(dim(i))}
  }
datasets <- list(dataset1, dataset2, dataset3, dataset4, dataset5, dataset6)
print_dimensions(datasets)


get_unique <- function(data_list) {
  for (i in seq_along(data_list)) {
    data_list[[i]] <- data_list[[i]] %>% distinct(Gene, .keep_all = TRUE)
  }
  return(data_list)
}
get_unique(datasets)
print_dimensions(datasets)

remove_duplicates <- function(data_list){
  for (i in data_list){
    i <- i[!duplicated(i[c('Gene')]), ]
    
  }
}

remove_duplicates(datasets)



make_numeric_columns <- function(dataset_list){
  for(i in dataset_list){
    i[,2:11] <- lapply(i[,2:11], as.numeric)}
}

make_numeric_columns(datasets)
str(dataset1)


make_row_names <- function(dataset_list){
  for (i in dataset_list){
    rownames(i) <- i$Gene}
}
make_row_names(datasets)
rownames(datasets[["dataset1"]])
################################### Still showing duplicates within the Gene column

remove_gene <- function(dataset_list) {
  for (i in seq_along(dataset_list)) {
    dataset_list[[i]] <- dataset_list[[i]] %>% filter(Gene != "GSN")
  }
  return(dataset_list)
}
################################### maybe i need to make it as a string or factor?

sapply(datasets, function(df) class(df$Gene))
gsn_present <- sapply(datasets, function(df) any(df$Gene == "GSN"))
print(gsn_present)
sum(sapply(datasets, function(df) sum(df$Gene == "GSN")))
datasets <- lapply(datasets, function(df) df[df$Gene != "GSN", ])
sum(sapply(datasets, function(df) sum(df$Gene == "GSN")))
print_dimensions(datasets)
################################ Removed! 
make_row_names(datasets)
datax <- dataset1[,5:10]
colnames(datax)
rownames(dataset1)
colnames(dataset1)



########## Data loaded again ########

read_and_set_rownames <- function(file) {
  df <- read.csv(file, header = TRUE)
    df <- as.data.frame(df)
    rownames(df) <- df$Gene
  return(df)
}
files <- c("cleaned1.csv", "cleaned2.csv", "cleaned3.csv", "cleaned4.csv", "cleaned5.csv", "cleaned6.csv")
datasets <- lapply(files, read_and_set_rownames)

dataset1 <- datasets[[1]]
dataset2 <- datasets[[2]]
dataset3 <- datasets[[3]]
dataset4 <- datasets[[4]]
dataset5 <- datasets[[5]]
dataset6 <- datasets[[6]]

################################### Heatmap ####################################
datasets <- list(dataset1, dataset2, dataset3, dataset4, dataset5, dataset6)

pipeline <- function(data_list){
  data_list_mean_centered <- list()
  for (i in data_list){
    i_subset <- i[,5:10]
    i_matrix <- as.matrix(i_subset)
    i_log <- log2(i_matrix + 1)
    i_mean_centered <- t(scale(t(i_log), scale = F))
    data_list_mean_centered[[i]] <- i_mean_centered
  }
  return(data_list_mean_centered)
}

processed_data <- pipeline(datasets)
dataset1_mean_centered <- processed_data[[1]]
head(dataset1_mean_centered)
dataset2_mean_centered <- processed_data[[2]]
dataset3_mean_centered <- processed_data[[3]]
dataset4_mean_centered <- processed_data[[4]]
dataset5_mean_centered <- processed_data[[5]]
dataset6_mean_centered <- processed_data[[6]]


#annotation_col <- data.frame(Group= c(sub("-.*", "", colnames(processed_data[[1]]))))
#annotation

data_vis <- list(dataset1_mean_centered, dataset2_mean_centered, dataset3_mean_centered, dataset4_mean_centered, dataset5_mean_centered, dataset6_mean_centered)

heatmap_maker <- function(data_list){
  for (i in data_list){
    i_annotation_col <- data.frame(Group= c(sub("-.*", "", colnames(i))))
    rownames(i_annotation_col) <- colnames(i)
    i_head <- pheatmap(i,
             cluster_cols = TRUE,  
             cluster_rows = TRUE,  
             scale = "row",  
             color = colorRampPalette(c("blue", "white", "red"))(100), 
             main = "Heatmap of Normalized Counts",  
             fontsize = 10,  
             fontsize_row = 8,  
             fontsize_col = 8,  
             #display_numbers = TRUE,  
             #number_color = "black",
             border_color = NA,  
             angle_col = 45,
             #rownames_force = TRUE, 
             annotation_col = i_annotation_col
    )
    ggsave(paste0("heatmap_", i, ".png"), i_head, width = 6, height = 4, units = "in")
  }
}

heatmap_maker(data_vis)

head(dataset1)
rownames(dataset1)
dataset1[,1] <- rownames(dataset1)
rownames(dataset1)
row.names(dataset1) <- dataset1[,1]
rownames(dataset1)
head(dataset1)
##################################### PCA ######################################
datasets <- list(dataset1, dataset2, dataset3, dataset4, dataset5, dataset6)

colnames(dataset1)
rownames(dataset1)
numberic_column_dataset1 <- dataset1[,5:10]
dataset1_numeric <- data.frame(lapply(numberic_column_dataset1, function(x) if(is.numeric(x)) as.integer(x) else x))

dataset1_numeric <- na.omit(dataset1_numeric)
colnames(dataset1_numeric)
#gr <- sub("\\..*", "", colnames(dataset1_numeric))
gr <- c(rep("CTRL", 3), rep("Eto", 3))
gr <- as.factor(gr)
gr
dataset1_numeric <- na.omit(dataset1_numeric)
colData <- data.frame(Group = gr)
rownames(colData) <- colnames(dataset1_numeric)
colnames(dataset1_numeric) <- gr
dim(dataset1_numeric)
cds <- DESeqDataSetFromMatrix(countData = dataset1_numeric, colData = colData, design = ~ Group)
cds <- DESeq(cds)
data.norm <- log2(1+counts(cds, normalized=F))
data.mean.center <- t(scale(t(data.norm), scale = F))
pc <- prcomp(data.mean.center)
pcr <- data.frame(pc$r)
pcr$group <- gr
pc_var <- pc$sdev^2 / sum(pc$sdev^2)

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



i <- dataset1
i <- i[, 5:10]
i <- na.omit(i)
i <- data.frame(lapply(i, function(x) if(is.numeric(x)) as.integer(x) else x))
i <- na.omit(i)
gr <- as.factor(gr)
i <- na.omit(i)
dds <- DESeqDataSetFromMatrix(countData = i,
                              colData = data.frame(Group=gr),
                              design= ~ Group)
dds <- DESeq(dds)
colData(dds)

data.norm <- log2(1+counts(dds, normalized=F))
data.mean.center <- t(scale(t(data.norm), scale = F))

pc <- prcomp(data.mean.center)
pcr <- data.frame(pc$r)
pcr$group <- gr
pc_var <- pc$sdev^2 / sum(pc$sdev^2)

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





datasets <- list(dataset1, dataset2, dataset3, dataset4, dataset5, dataset6)

pipeline_for_pca <- function(data_list){
  for (i in data_list){
    i <- i[, 5:10]
    i <- na.omit(i)
    i <- data.frame(lapply(i, function(x) if(is.numeric(x)) as.integer(x) else x))
    i <- na.omit(i)
    gr <- sub("\\..*", "", colnames(i))
    gr <- as.factor(gr)
    i <- na.omit(i)
    colData <- data.frame(Group = gr)
    rownames(colData) <- colnames(i)
    colnames(i) <- gr
    #cds <- DESeqDataSetFromMatrix(countData = i, colData = colData, design = ~ Group)
    cds <- DESeq(cds)
    data.norm <- log2(1+counts(cds, normalized=F))
    data.mean.center <- t(scale(t(data.norm), scale = F))
    pc <- prcomp(data.mean.center)
    pcr <- data.frame(pc$r)
    pcr$group <- gr
    pc_var <- pc$sdev^2 / sum(pc$sdev^2)
  }
}
pipeline_for_pca(datasets)












########## 

colnames(dataset1)

library(ggplot2)

pca_maker <- function(data_list) {
  for (i in seq_along(data_list)) {
    data <- data_list[[i]]
    data <- data[, 5:10]
    
    # Ensure all values are numeric and convert to integer, handle coercion issues
    data <- data.frame(lapply(data, function(x) {
      x <- as.numeric(as.character(x))
      if (any(is.na(x))) {
        stop("Non-numeric values detected in data.")
      }
      as.integer(x)
    }))
    
    # Remove rows with NA values
    data <- na.omit(data)
    
    # Remove columns with zero variance
    data <- data[, apply(data, 2, var, na.rm = TRUE) != 0]
    
    # Generate group labels from column names
    gr <- sub("\\..*", "", colnames(data))
    gr <- as.factor(gr)
    
    # Create column metadata
    colData <- data.frame(Group = gr)
    rownames(colData) <- colnames(data)
    
    # Log transform the data
    data.norm <- log2(1 + data)
    
    # Center the data by subtracting the mean
    data.mean.center <- t(scale(t(data.norm), scale = FALSE))
    
    # Perform PCA
    pc <- prcomp(data.mean.center)
    pcr <- data.frame(pc$x)
    pcr$group <- gr
    
    # Calculate variance explained by each PC
    pc_var <- pc$sdev^2 / sum(pc$sdev^2)
    
    # Plot PCA
    pca_plot <- ggplot(pcr, aes(PC1, PC2, color = group)) + 
      geom_point(size = 10, alpha = 0.6, shape = 19) + 
      theme_bw() + 
      theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
            axis.text = element_text(size = 12),
            axis.title = element_text(size = 14, face = "bold"),
            legend.title = element_text(size = 14, face = "bold"),
            legend.text = element_text(size = 12),
            legend.position = "right") + 
      labs(
        x = paste0("PC1 (", round(pc_var[1] * 100, 2), "% variance)"),
        y = paste0("PC2 (", round(pc_var[2] * 100, 2), "% variance)"),
        title = paste("PCA Plot - Dataset", i)
      )
    
    # Save the PCA plot
    ggsave(paste0("pca_", i, ".png"), pca_plot, width = 6, height = 4, units = "in")
  }
}

# Example usage with your datasets list
dataset <- list(dataset1, dataset2)
pca_maker(dataset)
