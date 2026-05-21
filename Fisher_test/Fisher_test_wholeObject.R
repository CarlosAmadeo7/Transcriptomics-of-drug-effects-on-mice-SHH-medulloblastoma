#-------------------------------------------------------------------------------
## Script to create input tables needed for Fisher enrichment test
# Loading Libraries
suppressPackageStartupMessages({
  ;library(ggplot2);library(ggrepel);library(ggpubr);library(data.table);library(RColorBrewer)
  ;library(tidyverse);library(preprocessCore);library(future.apply);library(DESeq2);library(pheatmap)
  ;library(sva);library(viridis);library(limma);library(emmeans);library(broom);library(janitor)
  ;library(tidyplots);library(dplyr);library(writexl);library(scToppR);library(fgsea)
  ;library(clusterProfiler);library(org.Hs.eg.db);library(dplyr);library(openxlsx);library(readxl);library(Seurat)})

# Set directory of work ( files in here)
#setwd("MUSC/Bioinformatics/BMS.Ribo.Combo/")

## Loading Single cell Object ( previosuly analyzed and annotated)
load("C:/Users/amade/Documents/Pritika's_project_BMS_Ribo_Combo/")
DimPlot(seuObject_slim_nodoub, label = T)

# Creating a generalized annotation for Fisher_test
Idents(seuObject_slim_nodoub)<- seuObject_slim_nodoub$SCT_snn_res.0.8
unique(Idents(seuObject_slim_nodoub))
DefaultAssay(seuObject_slim_nodoub)<- "RNA"

cluster_annotations <- c(
  "0"  = "Differentiating tumor cells",
  "1"  = "Cycling tumor cells",
  "2"  = "Differentiating tumor cells",
  "3"  = "GNP-like cells",
  "4"  = "GNP-like cells",
  "5"  = "GNP-like cells",
  "6"  = "Cycling tumor cells",
  "7"  = "Cycling tumor cells",
  "8"  = "Differentiating tumor cells",
  "9"  = "Cycling tumor cells",
  "10" = "Oligodendrocyte-precursor cells",
  "11" = "Mature neurons",
  "12" = "Pvalb+ tumor cells",
  "13" = "Astrocytes",
  "14" = "Myeloid cells",
  "15" = "Fibroblasts",
  "16" = "Endothelial Cells"
)
Idents(seuObject_slim_nodoub)<-seuObject_slim_nodoub$SCT_snn_res.0.8
seuObject_slim_nodoub@meta.data$celltype1<- cluster_annotations[as.character(Idents(seuObject_slim_nodoub))]
DimPlot(seuObject_slim_nodoub, group.by = "celltype1")

# Function to find DEGs/ cluster
DefaultAssay(seuObject_slim_nodoub)<- "RNA"
Idents(seuObject_slim_nodoub)<- seuObject_slim_nodoub$celltype1
DEGs_Object <- Reduce("rbind",lapply(unique(seuObject_slim_nodoub$celltype1), function(x) {
  Markers <- FindMarkers(seuObject_slim_nodoub, ident.1 = x, ident.2 = NULL, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.5, pseudocount.use = 0.1)
  Markers <- Markers[which(Markers$p_val_adj < 0.05),]
  Markers$gene <- rownames(Markers)
  Markers$cluster<- rep(paste(x), nrow(Markers))
  return(Markers)}))

# Create GeneSets
gene.clusters<- DEGs_Object %>% dplyr::select(gene, cluster)
colnames(gene.clusters)= c("Gene","Clusters")
GeneSets <- split(gene.clusters, gene.clusters$Clusters)
save(GeneSets, file= "utils/geneset/GeneSets_Palbociclib_DMSO.RData") # Save

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
RData <-load("Enrichments/STATS/Pritika_project_markers_Ribo_Stats.RData")
print(FisherMat$Astrocytes)
#significant_clusters <- lapply(FisherMat, function(x) x[order(x[, "P_val"]), ])
for (name in names(FisherMat)) {cat("\nGene set:", name, "\n")print(FisherMat[[name]])}
  