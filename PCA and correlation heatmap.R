setwd("/Users/mortezaabyadeh/Desktop")
library(data.table)
library(ggplot2)
library(DESeq2)
library(openxlsx) #don't need here but in case for other analysis I kept it here
library(readxl)
library(reshape2)
data <- read_excel("B_file.xlsx")
dim(data)

class(data)
colnames(data)
data[,1] <- rownames(data)
data <- data[,-1]
gr <- c("B_IDO-F-B", "B_IDO-F-B-ctrl", "B_IDO-M-B", "B_IDO-M-B-ctrl", "B_IDO1-F-B", "B_IDO1-F-B-ctrl", "B_IDO1-M-B", "B_IDO1-M-B-ctrl", "B_IDO2-F-B",
        "B_IDO2-F-B-ctrl", "B_IDO2-M-B", "B_IDO2-M-B-ctrl", "U_IDO-F-B", "U_IDO-F-B-ctrl", "U_IDO-M-B", "U_IDO-M-B-ctrl", "U_IDO1-F-B", "U_IDO1-F-B-ctrl", "U_IDO1-M-B", 
        "U_IDO1-M-B-ctrl", "U_IDO2-F-B", "U_IDO2-F-B-ctrl", "U_IDO2-M-B", "U_IDO2-M-B-ctrl")

gr <- as.factor(gr)
data <- na.omit(data)

# Convert list to data frame
data <- as.data.frame(data)

# Check dimensions after conversion
dim(data)

# Proceed with selecting numeric columns
numeric_columns <- sapply(data, is.numeric)

data_numeric <- data[, which(numeric_columns)]
data_numeric <- round(data_numeric)
colData <- data.frame(Group = gr)
rownames(colData) <- colnames(data)

summary(data_numeric)

cds <- DESeqDataSetFromMatrix(countData = data_numeric, colData = colData, design = ~ Group)
cds <- DESeq(cds)
head(data)

cds <- estimateSizeFactors(cds)
data.norm <- log2(1+counts(cds, normalized=T))
head(data.norm)
boxplot(data.norm)

data.mean.center <- t(scale(t(data.norm), scale = F))


head(data.mean.center)

pc <- prcomp(data.mean.center)
plot(pc)
head(pc$rotation)
pcr <- data.frame(pc$r)
head(pcr)
head(data)
pcr$group <- gr
head(pcr)


head(data)

# write.xlsx(data, "Data_R_PCA.xlsx")


#################################### PCA ######################################


ggplot(pcr, aes(PC1, PC2, color = group, shape = group)) + 
  geom_point(size = 5, alpha = 0.9) + 
  scale_shape_manual(values = 1:24) +  # Manually specify shapes if needed
  theme_bw() + 
  theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.position = "right"
  )


ggplot(pcr, aes(PC1, PC2, color = group)) + 
  geom_point(aes(shape = "Group"), size = 5, alpha = 0.9) + 
  scale_shape_manual(values = c("Group" = 19)) +  
  theme_bw() + 
  theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.position = "right"
  ) + guides(color = guide_legend(override.aes = list(shape = c(19))))


##################################### Correlation heat map 

cor_matrix <- cor(data.mean.center)
cor_melted <- melt(cor_matrix)




ggplot(cor_melted, aes(Var2, Var1, fill = value, label = round(value, 2))) +
  geom_tile(color = "white") +
  geom_text(size = 3, color = "black") +  
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                       limits = c(-1, 1), name = "Correlation") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    axis.title = element_blank()
  ) +
  labs(
    title = "Correlation Heatmap",
    x = "Variables",
    y = "Variables"
  )


######################## Upper correlation heat map
cor_upper <- cor_matrix
cor_upper[lower.tri(cor_upper, diag = TRUE)] <- NA  


cor_upper_melted <- melt(cor_upper, na.rm = TRUE)


ggplot(cor_upper_melted, aes(Var2, Var1, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                       limits = c(-1, 1), name = "Correlation") +
  geom_text(aes(label = round(value, 2)), color = "black", size = 3) +  
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    axis.title = element_blank()
  ) +
  labs(
    title = "Upper Triangle Correlation Heatmap",
    x = "Variables",
    y = "Variables"
  ) +
  coord_fixed()  





############################### For Liver Dataset ############################## 

data <- read_excel("L_file.xlsx")
dim(data)

class(data)
colnames(data)
data[,1] <- rownames(data)
data <- data[,-1]
gr <- c("B_IDO-F-L", "B_IDO-F-L-ctrl", "B_IDO-M-L", "B_IDO-M-L-ctrl", "B_IDO1-F-L", "B_IDO1-F-L-ctrl", "B_IDO1-M-L", "B_IDO1-M-L-ctrl", "B_IDO2-F-L",
        "B_IDO2-F-L-ctrl", "B_IDO2-M-L", "B_IDO2-M-L-ctrl", "U_IDO-F-L", "U_IDO-F-L-ctrl", "U_IDO-M-L", "U_IDO-M-L-ctrl", "U_IDO1-F-L", "U_IDO1-F-L-ctrl", "U_IDO1-M-L", 
        "U_IDO1-M-L-ctrl", "U_IDO2-F-L", "U_IDO2-F-L-ctrl", "U_IDO2-M-L", "U_IDO2-M-L-ctrl")

gr <- as.factor(gr)
data <- na.omit(data)

# Convert list to data frame
data <- as.data.frame(data)

# Check dimensions after conversion
dim(data)

# Proceed with selecting numeric columns
numeric_columns <- sapply(data, is.numeric)

data_numeric <- data[, which(numeric_columns)]
data_numeric <- round(data_numeric)
colData <- data.frame(Group = gr)
rownames(colData) <- colnames(data)

summary(data_numeric)

cds <- DESeqDataSetFromMatrix(countData = data_numeric, colData = colData, design = ~ Group)
cds <- DESeq(cds)
head(data)

cds <- estimateSizeFactors(cds)
data.norm <- log2(1+counts(cds, normalized=T))
head(data.norm)
boxplot(data.norm)

data.mean.center <- t(scale(t(data.norm), scale = F))

pc <- prcomp(data.mean.center)
plot(pc)
head(pc$rotation)
pcr <- data.frame(pc$r)
head(pcr)
head(data)
pcr$group <- gr
head(pcr)


head(data)

# write.xlsx(data, "Data_R_PCA.xlsx")



ggplot(pcr, aes(PC1, PC2, color = group, shape = group)) + 
  geom_point(size = 5, alpha = 0.9) + 
  scale_shape_manual(values = 1:24) +  # Manually specify shapes if needed
  theme_bw() + 
  theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.position = "right"
  )


ggplot(pcr, aes(PC1, PC2, color = group)) + 
  geom_point(aes(shape = "Group"), size = 5, alpha = 0.9) + 
  scale_shape_manual(values = c("Group" = 19)) +  
  theme_bw() + 
  theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.position = "right"
  ) + guides(color = guide_legend(override.aes = list(shape = c(19))))



####################################### Correlation heat map ####################
cor_matrix <- cor(data.mean.center)
cor_melted <- melt(cor_matrix)




ggplot(cor_melted, aes(Var2, Var1, fill = value, label = round(value, 2))) +
  geom_tile(color = "white") +
  geom_text(size = 3, color = "black") +  
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                       limits = c(-1, 1), name = "Correlation") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    axis.title = element_blank()
  ) +
  labs(
    title = "Correlation Heatmap",
    x = "Variables",
    y = "Variables"
  )

########################################## Lowe side
cor_lower <- cor_matrix
cor_lower[upper.tri(cor_lower, diag = TRUE)] <- NA  


cor_lower_melted <- melt(cor_lower, na.rm = TRUE)

ggplot(cor_lower_melted, aes(Var2, Var1, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                       limits = c(-1, 1), name = "Correlation") +
  geom_text(aes(label = round(value, 2)), color = "black", size = 3) +  
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    axis.title = element_blank()
  ) +
  labs(
    title = "Triangular Correlation Heatmap",
    x = "Variables",
    y = "Variables"
  ) +
  coord_fixed()  



################################# Upper side
cor_upper <- cor_matrix
cor_upper[lower.tri(cor_upper, diag = TRUE)] <- NA  


cor_upper_melted <- melt(cor_upper, na.rm = TRUE)


ggplot(cor_upper_melted, aes(Var2, Var1, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                       limits = c(-1, 1), name = "Correlation") +
  geom_text(aes(label = round(value, 2)), color = "black", size = 3) +  
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    axis.title = element_blank()
  ) +
  labs(
    title = "Upper Triangle Correlation Heatmap",
    x = "Variables",
    y = "Variables"
  ) +
  coord_fixed()  

