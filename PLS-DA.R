
library(DESeq2)
library(readxl)
library(dplyr)
library(pheatmap)
library(data.table)
library(ggplot2)
library(openxlsx) # not needed here but kept for other analysis
library(plotly)
library(htmlwidgets)
library(RColorBrewer) 
library(ggrepel)
library(tidyverse)
library(reshape2)
library(shiny)
library(DOSE)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
#library(pathfindR)
library(enrichR)
library(biomaRt)

setwd("/Users/mortezaabyadeh/Desktop/senescence omics data/used for app")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("mixOmics")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ropls")


library("ropls")
library(mixOmics)
library(ggplot2)

transcriptome_Eto_vs_CTRL <- read_excel("clean_transcriptome_Eto_vs_CTRL.xlsx")
transcriptome_FTY_vs_CTRL <- read_excel("clean_transcriptome_FTY_vs_CTRL.xlsx")
transcriptome_FTY_Eto_vs_CTRL <- read_excel("clean_transcriptome_FTY_Eto_vs_CTRL.xlsx")
transcriptome_FTY_vs_Eto <- read_excel("clean_transcriptome_FTY_vs_Eto.xlsx")
transcriptome_FTY_Eto_vs_Eto <- read_excel("clean_transcriptome_FTY_Eto_vs_Eto.xlsx")
transcriptome_FTY_Eto_vs_FTY <- read_excel("clean_transcriptome_FTY_Eto_vs_FTY.xlsx")



proteome_Eto_vs_CTRL <- read_excel("set1.xlsx")
proteome_FTY_vs_CTRL<- read_excel("set2.xlsx")
proteome_FTY_Eto_vs_CTRL <- read_excel("set3.xlsx")
proteome_FTY_vs_Eto <- read_excel("set4.xlsx")
proteome_FTY_Eto_vs_Eto <- read_excel("set5.xlsx")
proteome_FTY_Eto_vs_FTY <- read_excel("set6.xlsx")


colnames(proteome_Eto_vs_CTRL)[which(names(proteome_Eto_vs_CTRL) == "p-adj")] <- "p_adj"
colnames(proteome_FTY_vs_CTRL)[which(names(proteome_FTY_vs_CTRL) == "p-adj")] <- "p_adj"
colnames(proteome_FTY_Eto_vs_CTRL)[which(names(proteome_FTY_Eto_vs_CTRL) == "p-adj")] <- "p_adj"
colnames(proteome_FTY_vs_Eto)[which(names(proteome_FTY_vs_Eto) == "p-adj")] <- "p_adj"
colnames(proteome_FTY_Eto_vs_Eto)[which(names(proteome_FTY_Eto_vs_Eto) == "p-adj")] <- "p_adj"
colnames(proteome_FTY_Eto_vs_FTY)[which(names(proteome_FTY_Eto_vs_FTY) == "p-adj")] <- "p_adj"

print(round(-log10(min(proteome_Eto_vs_CTRL$p_adj))))
min(proteome_Eto_vs_CTRL$p_adj)


datasets <- list(
  transcriptome = list(
    Eto_vs_CTRL = transcriptome_Eto_vs_CTRL,
    FTY_vs_CTRL = transcriptome_FTY_vs_CTRL,
    FTY_Eto_vs_CTRL = transcriptome_FTY_Eto_vs_CTRL,
    FTY_vs_Eto = transcriptome_FTY_vs_Eto,
    FTY_Eto_vs_Eto = transcriptome_FTY_Eto_vs_Eto,
    FTY_Eto_vs_FTY = transcriptome_FTY_Eto_vs_FTY
  ),
  proteome = list(
    Eto_vs_CTRL = proteome_Eto_vs_CTRL,
    FTY_vs_CTRL = proteome_FTY_vs_CTRL,
    FTY_Eto_vs_CTRL = proteome_FTY_Eto_vs_CTRL,
    FTY_vs_Eto = proteome_FTY_vs_Eto,
    FTY_Eto_vs_Eto = proteome_FTY_Eto_vs_Eto,
    FTY_Eto_vs_FTY = proteome_FTY_Eto_vs_FTY
  )
)


plsda_maker <- function(df){
  df_filtered <- proteome_Eto_vs_CTRL[proteome_Eto_vs_CTRL$`p_adj` < 0.05, ]
  numeric_columns <- df_filtered[,2:7]
  head(df_filtered)
  dataset_numeric <- data.frame(lapply(numeric_columns, function(x) if(is.numeric(x)) as.integer(x) else x))
  dim(dataset_numeric)
  dataset_numeric <- na.omit(dataset_numeric)
  dim(dataset_numeric)
  head(dataset_numeric)
  group <- sub("\\..*", "", colnames(dataset_numeric))
  group <- as.factor(group)
  group
  levels = unique(group)
  group
  colnames(dataset_numeric) <- unlist(lapply(levels, function(g) {
    paste0(g, "-", seq(sum(group == g)))
  }))
  
  colnames(dataset_numeric)
  
  plsda_model <- plsda(t(dataset_numeric), group, ncomp = 2)
  
  length(group)
  nrow(plsda_model$X)
  str(plsda_model)
  
  plotIndiv(plsda_model, comp = c(1, 2), 
            group = group, 
            ind.names = FALSE, 
            ellipse = FALSE, 
            title = "PLS-DA Plot")
}

plsda_maker(proteome_Eto_vs_CTRL)



