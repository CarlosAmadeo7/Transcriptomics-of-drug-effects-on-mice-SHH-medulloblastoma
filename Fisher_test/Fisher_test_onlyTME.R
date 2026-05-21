#-------------------------------------------------------------------------------
## Script to create input tables needed for Fisher enrichment test (from TME)
# Loading Libraries
suppressPackageStartupMessages({
  ;library(ggplot2);library(ggrepel);library(ggpubr);library(data.table);library(RColorBrewer)
  ;library(tidyverse);library(preprocessCore);library(future.apply);library(DESeq2);library(pheatmap)
  ;library(sva);library(viridis);library(limma);library(emmeans);library(broom);library(janitor)
  ;library(tidyplots);library(dplyr);library(writexl);library(scToppR);library(fgsea)
  ;library(clusterProfiler);library(org.Hs.eg.db);library(dplyr);library(openxlsx);library(readxl);library(Seurat)})

# Set directory of work ( files in here)
#setwd("MUSC/Bioinformatics/BMS.Ribo.Combo/")

## Loading Single cell Object 
#load("C:/Users/amade/Documents/Pritika's_project_BMS_Ribo_Combo/  ")
DimPlot(TME, label = T)
DimPlot(TME, group.by = "celltype2")

Idents(TME)<- TME$celltype2
DefaultAssay(TME)<- "RNA"

# Function to find DEGs/ cluster
DefaultAssay(TME)<- "RNA"
Idents(TME)<- TME$celltype2
DEGs_Object <- Reduce("rbind",lapply(unique(TME$celltype2), function(x) {
  Markers <- FindMarkers(TME, ident.1 = x, ident.2 = NULL, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.5, pseudocount.use = 0.1)
  Markers <- Markers[which(Markers$p_val_adj < 0.05),]
  Markers$gene <- rownames(Markers)
  Markers$cluster<- rep(paste(x), nrow(Markers))
  return(Markers)}))

# Create GeneSets
gene.clusters<- DEGs_Object %>% dplyr::select(gene, cluster)
colnames(gene.clusters)= c("Gene","Clusters")
GeneSets <- split(gene.clusters, gene.clusters$Clusters)
save(GeneSets, file= "utils/geneset/GeneSets_TME_Palbociclib_DMSO.RData")

###------------------
## Loading the DEG objects from Combo, Ribo , BMS vs DMSO ( bulk RNA-seq)
Combo<-read_xlsx("C:/Users/amade/Documents/Pritika's_project_BMS_Ribo_Combo/edgeR_new_exploration/deg/Combo_DMSO/CTX_DGE_shrink.xlsx")
BMS<-read_xlsx("C:/Users/amade/Documents/Pritika's_project_BMS_Ribo_Combo/edgeR_new_exploration/deg/BMS_DMSO/CTX_DGE_shrink.xlsx")
Ribo<- read_xlsx("C:/Users/amade/Documents/Pritika's_project_BMS_Ribo_Combo/edgeR_new_exploration/deg/Ribo_DMSO/CTX_DGE_shrink.xlsx")

# Temporary frames for each comparison and save as txt files
# Combo vs DMSO 
df <- Combo %>%
  mutate(LOG=-log10(padj), ABS=abs(log2FoldChange))%>%
  mutate(Threshold = if_else(padj<0.05 & ABS>0.3, "TRUE", "FALSE"))%>%
  mutate(Direction = case_when(log2FoldChange>0.3 & padj<0.05 ~ "Up", log2FoldChange< -0.3 & padj<0.05 ~ "Down")) %>%
  dplyr::select(Gene, Direction)
all <- data.frame(Gene=df$Gene,Direction="All")
deg = rbind(df, all)
write.table(deg, "deg.Combo.DMSO.txt", sep="\t", quote = F)

# BMS vs DMSO 
df <- BMS %>%
  mutate(LOG=-log10(padj), ABS=abs(log2FoldChange))%>%
  mutate(Threshold = if_else(padj<0.05 & ABS>0.3, "TRUE", "FALSE"))%>%
  mutate(Direction = case_when(log2FoldChange>0.3 & padj<0.05 ~ "Up", log2FoldChange< -0.3 & padj<0.05 ~ "Down")) %>%
  dplyr::select(Gene, Direction)
all <- data.frame(Gene=df$Gene,Direction="All")
deg = rbind(df, all)
write.table(deg, "deg.BMS.DMSO.txt", sep="\t", quote = F)

# Ribo vs DMSO 
df <- Ribo %>%
  mutate(LOG=-log10(padj), ABS=abs(log2FoldChange))%>%
  mutate(Threshold = if_else(padj<0.05 & ABS>0.3, "TRUE", "FALSE"))%>%
  mutate(Direction = case_when(log2FoldChange>0.3 & padj<0.05 ~ "Up", log2FoldChange< -0.3 & padj<0.05 ~ "Down")) %>%
  dplyr::select(Gene, Direction)
all <- data.frame(Gene=df$Gene,Direction="All")
deg = rbind(df, all)
write.table(deg, "deg.Ribo.DMSO.txt", sep="\t", quote = F)

###----------------------------------------------------------
## Opening the RDATA sets to inspect
RData<-load("sC.comparasion/enrichments/STATS/Medullo_Markers_Stats.RData")
head(RData)
View(FisherMat)
print(FisherMat$`Cluster 1 `)
significant_clusters <- lapply(FisherMat, function(x) x[order(x[, "P_val"]), ])
