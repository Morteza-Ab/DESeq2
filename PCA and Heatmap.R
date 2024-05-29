setwd("/Users/mortezaabyadeh/Desktop")
library(data.table)
library(ggplot2)
library(DESeq2)
library(openxlsx) #don't need here but in case for other analysis I kept it here
library(readxl)
data <- read_excel("combined_proteome_transcriptome_data.xlsx")
class(data)
colnames(data)

data[,1] <- rownames(data)
data <- data[,-1]
gr <- c(rep("proteome-lsm5", 3), rep("proteome-ctrl", 3), rep("transcriptome-lsm5", 3), rep("transcriptome-ctrl", 3))
gr <- as.factor(gr)

data <- na.omit(data)

######## ######. You can use the following code to convert to integer
data <- data[, sapply(data, is.numeric)] <- lapply(data[, sapply(data, is.numeric)], as.integer)
dim(data)
head(data, 5)

####### ###### I got some error making NA while trying to convert floating to integer using the above code
#### So the below code can work better in such cases

cleaned_data <- as.data.frame(lapply(cleaned_data, function(x) {
  x <- as.numeric(as.character(x))
  x[is.na(x)] <- 0  # Replace NAs with 0
  x[x > .Machine$integer.max] <- .Machine$integer.max  # Cap values at integer max
  x[x < .Machine$integer.min] <- .Machine$integer.min  # Cap values at integer min
  as.integer(x)}))

################################## I had NA in the data: if (any(is.na(data))) {
# print("There are NA values in the count matrix.")

total_na <- sum(is.na(cleaned_data))
print(paste("Total number of NA values:", total_na)) 

cleaned_data <- na.omit(data)
dim(cleaned_data)
write.xlsx(cleaned_data, "cleaned1_combined_data.xlsx")
data <- na.omit(data)


# Ensure that the number of columns and rows match
dim(cleaned_data)  # Expecting (genes x samples)
dim(colData)       # Expecting (samples x metadata_columns)


# Check sample names
colnames(cleaned_data)
rownames(colData)

# Subset colData to match the samples in cleaned_data if necessary
## colData <- colData[colnames(cleaned_data), , drop = FALSE]



colData <- data.frame(Group = gr)
rownames(colData) <- colnames(cleaned_data)
cds <- DESeqDataSetFromMatrix(countData = cleaned_data, colData = colData, design = ~ Group)
cds <- DESeq(cds)
head(data)


data.norm <- log2(1+counts(cds, normalized=T))
head(data.norm)
boxplot(data.norm)

###Function_I made that_can help to more discriminate between data in pca plot-You may do not like it!
data.mean.center <- t(scale(t(data.norm), scale = F))


###pcr
pc <- prcomp(data.mean.center)
plot(pc)
head(pc$rotation)
pcr <- data.frame(pc$r)
head(pcr)
head(data)
pcr$group <- gr
head(pcr)


library(ggplot2)

ggplot(pcr, aes(PC1, PC2, color = group, shape = group)) + 
  geom_point(size = 5, alpha = 0.9) + 
  theme_bw() + 
  theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.position = "right"
  )

ggplot(pcr, aes(PC1, PC2, color=group)) + geom_point(size=5, alpha=0.9, shape = 10) + theme_bw() + theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
                                                                                                         axis.text = element_text(size = 12),
                                                                                                         axis.title = element_text(size = 14, face = "bold"),
                                                                                                         legend.title = element_text(size = 14, face = "bold"),
                                                                                                         legend.text = element_text(size = 12),
                                                                                                         legend.position = "right")
cbPalette <- c("#999999", "#E69F00", "#56B4E9")

custom_palette <- c("lightgreen", "gold", "lightblue","purple")

ggplot(pcr, aes(PC1, PC2, color = group, shape = group)) +
  geom_point(size = 9, alpha = 0.9, position = position_jitter(width = 0.1, height = 0.1)) +  
  scale_color_manual(values = custom_palette) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15, face = "bold"),
        legend.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 15),
        legend.position = "right")

################################################## Heatmap #############################################
library(gplots)

normalized_counts <- counts(cds, normalized = TRUE)

heatmap(normalized_counts, 
        Colv = NA,  # Turn off column clustering
        Rowv = NA,  # Turn off row clustering
        scale = "row",  # Scale rows (genes)
        col = colorRampPalette(c("blue", "white", "red"))(100),  # Color palette
        main = "Heatmap of Normalized Counts",  # Title
        xlab = "Samples",  # X-axis label
        ylab = "Genes"  # Y-axis label
)


library(pheatmap)
pheatmap(data.mean.center)
?pheatmap
pheatmap(data.norm)
pheatmap(cor(data.mean.center))
pheatmap(cor(data.mean.center))
head(data.mean.center)

####################making heatmap only for proteome data
dim(data.norm)
data.norm1 <- data.norm[, 1:6]
dim(data.norm1)
head(data.norm1)
pheatmap(data.norm1)
pheatmap(data.mean.center1)
data.mean.center1 <- t(scale(t(data.norm1), scale = F))
pheatmap(cor(data.norm1))
pheatmap(cor(data.mean.center1))


library(gplots)

# normalized_counts1 <- counts(data.mean.center1, normalized = TRUE)

heatmap(data.mean.center1, 
        Colv = NA,  # Turn off column clustering
        Rowv = NA,  # Turn off row clustering
        scale = "row",  # Scale rows (genes)
        col = colorRampPalette(c("blue", "white", "red"))(100),  # Color palette
        main = "Heatmap of Normalized Counts",  # Title
        xlab = "Samples",  # X-axis label
        ylab = "Genes"  # Y-axis label
)



### Pheatmap
library(pheatmap)
data.mean.center1 <- matrix(rnorm(100), nrow = 10, ncol = 10)
rownames(data.mean.center1) <- paste("Gene", 1:10, sep = "")
colnames(data.mean.center1) <- paste("Sample", 1:10, sep = "")

# Create heatmap using pheatmap and correlation heatmap #####
pheatmap(data.mean.center1,
         cluster_cols = FALSE,  # Turn off column clustering
         cluster_rows = FALSE,  # Turn off row clustering
         scale = "row",  # Scale rows (genes)
         color = colorRampPalette(c("blue", "white", "red"))(100),  # Color palette
         main = "Heatmap of Normalized Counts",  # Title
         fontsize = 10,  # Font size for general text
         fontsize_row = 8,  # Font size for row labels
         fontsize_col = 8,  # Font size for column labels
         display_numbers = TRUE,  # Display cell values
         number_color = "black",  # Color for numbers
         border_color = NA,  # No border around cells
         angle_col = 45,  # Rotate column labels for readability
         annotation_col = data.frame(Group = rep(c("A", "B"), each = 5))  # Add sample annotations
)


data.mean.center1 <- matrix(rnorm(100), nrow = 10, ncol = 10)
rownames(data.mean.center1) <- paste("Gene", 1:10, sep = "")
colnames(data.mean.center1) <- paste("Sample", 1:10, sep = "")

# Create an annotation data frame
annotation_col <- data.frame(Group = rep(c("A", "B"), each = 5))
rownames(annotation_col) <- colnames(data.mean.center1)

# Create heatmap using pheatmap
pheatmap(data.mean.center1,
         cluster_cols = FALSE,  # Turn off column clustering
         cluster_rows = FALSE,  # Turn off row clustering
         scale = "row",  # Scale rows (genes)
         color = colorRampPalette(c("blue", "white", "red"))(100),  # Color palette
         main = "Heatmap of Normalized Counts",  # Title
         fontsize = 10,  # Font size for general text
         fontsize_row = 8,  # Font size for row labels
         fontsize_col = 8,  # Font size for column labels
         display_numbers = TRUE,  # Display cell values
         number_color = "black",  # Color for numbers
         border_color = NA,  # No border around cells
         angle_col = 45,  # Rotate column labels for readability
         annotation_col = annotation_col  # Add sample annotations
)

################## Heatmap for only proteome data ################

data <- read_excel("combined_proteome_transcriptome_data.xlsx")
data_proteome <- data[,1:7]
head(data_proteome)
dim(data_proteome)

library(pheatmap)
rownames(data_proteome) <- data_proteome$`Gene name`
data_matrix <- as.matrix(data_proteome[, -1])

# Create sample annotations (assuming "lsm5" and "ctrl" groups)
annotation_col <- data.frame(Group = c("lsm5", "lsm5", "lsm5", "ctrl", "ctrl", "ctrl"))
rownames(annotation_col) <- colnames(data_matrix)

# Create heatmap using pheatmap
pheatmap(data_matrix,
         cluster_cols = FALSE,  # Turn off column clustering
         cluster_rows = FALSE,  # Turn off row clustering
         scale = "row",  # Scale rows (genes)
         color = colorRampPalette(c("blue", "white", "red"))(100),  # Color palette
         main = "Heatmap of Normalized Counts",  # Title
         fontsize = 10,  # Font size for general text
         fontsize_row = 8,  # Font size for row labels
         fontsize_col = 8,  # Font size for column labels
         #display_numbers = TRUE,  # Display cell values
         #number_color = "black",  # Color for numbers
         border_color = NA,  # No border around cells
         angle_col = 45,  # Rotate column labels for readability
         annotation_col = annotation_col  # Add sample annotations
)

###### adding clustering to heatmap

pheatmap(data_matrix,
         cluster_cols = TRUE,  # Enable column clustering
         cluster_rows = TRUE,  # Enable row clustering
         scale = "row",  # Scale rows (genes)
         color = colorRampPalette(c("blue", "white", "red"))(100),  # Color palette
         main = "Heatmap of Normalized Counts",  # Title
         fontsize = 10,  # Font size for general text
         fontsize_row = 8,  # Font size for row labels
         fontsize_col = 8,  # Font size for column labels
         # display_numbers = TRUE,  # Display cell values (commented out)
         # number_color = "black",  # Color for numbers (commented out)
         border_color = NA,  # No border around cells
         angle_col = 45,  # Rotate column labels for readability
         annotation_col = annotation_col  # Add sample annotations
)
  
###### adding gene name

pheatmap(data_matrix,
         cluster_cols = TRUE,  # Enable column clustering
         cluster_rows = TRUE,  # Enable row clustering
         scale = "row",  # Scale rows (genes)
         color = colorRampPalette(c("blue", "white", "red"))(100),  # Color palette
         main = "Heatmap of Normalized Counts",  # Title
         fontsize = 10,  # Font size for general text
         fontsize_row = 20,  # Font size for row labels
         fontsize_col = 12,  # Font size for column labels
         # display_numbers = TRUE,  # Display cell values (commented out)
         # number_color = "black",  # Color for numbers (commented out)
         border_color = NA,  # No border around cells
         angle_col = 45,  # Rotate column labels for readability
         annotation_col = annotation_col,  # Add sample annotations
         rownames_force = TRUE  # Force display of row names
)


##################### make only for first 10 genes

data <- read_excel("combined_proteome_transcriptome_data.xlsx")
data_proteome <- data[ 1:11, 1:7]
head(data_proteome)
dim(data_proteome)
rownames(data_proteome) <- data_proteome$`Gene name`
rownames(data_proteome)
data_matrix <- as.matrix(data_proteome[, -1])

head(data_proteome)

# Create sample annotations (assuming "lsm5" and "ctrl" groups)
annotation_col <- data.frame(Group = c("lsm5", "lsm5", "lsm5", "ctrl", "ctrl", "ctrl"))
rownames(annotation_col) <- colnames(data_matrix)
rownames(data_matrix) <- rownames(data_proteome)


pheatmap(data_matrix,
         cluster_cols = TRUE,  # Enable column clustering
         cluster_rows = TRUE,  # Enable row clustering
         scale = "row",  # Scale rows (genes)
         color = colorRampPalette(c("blue", "white", "red"))(100),  # Color palette
         main = "Heatmap of Normalized Counts",  # Title
         fontsize = 10,  # Font size for general text
         fontsize_row = 20,  # Font size for row labels
         fontsize_col = 12,  # Font size for column labels
         # display_numbers = TRUE,  # Display cell values (commented out)
         # number_color = "black",  # Color for numbers (commented out)
         border_color = NA,  # No border around cells
         angle_col = 45,  # Rotate column labels for readability
         annotation_col = annotation_col,  # Add sample annotations
         rownames_force = TRUE  # Force display of row names
)

