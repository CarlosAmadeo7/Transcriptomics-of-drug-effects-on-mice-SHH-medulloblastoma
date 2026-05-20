#' ---
#' title: "Vehicle sample analysis for bulk RNA seq Fisher-test enrichment. Data obtained from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE150579"
#' author: "Carlos Alfaro"
#' date: "2025-04-09"
#' output: html_document
#' ---
#' 
## ----setup, include=FALSE------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

#' 
## ------------------------------------------------------------------------------------------------------------------
library(SingleCellExperiment);library(Seurat);library(tidyverse);library(Matrix)
;library(scales);library(cowplot);library(RCurl); library(pheatmap)
;library(openxlsx);library(knitr);library(SeuratWrappers);library(scuttle); 
library(presto); library(RColorBrewer); library(sessioninfo)

#' 
#' #---------------
#' # Pre-processing
## ----pre-processing------------------------------------------------------------------------------------------------
# Loading original Seurat object from GSE150579
object<- readRDS("data/GSE150579_GSMO_Working.rds/GSE150579_GSMO_Working.rds")
# Update Seurat Object
set.seed(7081998)
GT <- UpdateSeuratObject(object = object)

#' 
## ------------------------------------------------------------------------------------------------------------------
# Adding metadata
metadata_mapping <- data.frame(
  SampleID = c( "G6" ,   "G7",    "G8",    "G9",    "G10",   "G11"),
  Treatment = c("Vehicle", "Vehicle", "Vehicle", "Vehicle", "Vehicle", "Vehicle"))
# Sample name modification 
cell_sample_ids <- sapply(strsplit(colnames(GT), "_"), `[`, 1)
cell_metadata <- data.frame(Treatment = metadata_mapping$Treatment[match(cell_sample_ids, metadata_mapping$SampleID)])
rownames(cell_metadata) <- colnames(GT)
# Merge metadata to original object 
GT <- AddMetaData(GT, metadata = cell_metadata)

#' 
## ----mitochondrial, ribosomal percentages--------------------------------------------------------------------------
#percent.mt
GT[['pMito']] <- PercentageFeatureSet(GT, pattern = '^mt-')
if(sum(GT[['pMito']]) == 0) {
  GT[['pMito']] <- PercentageFeatureSet(GT, pattern = '^MT-')}
if(sum(GT[['pMito']]) == 0) {
  GT[['pMito']] <- PercentageFeatureSet(GT, pattern = '^Mt-')}
# percent.ribo
GT[["pRibo"]] <- PercentageFeatureSet(GT, pattern = "^Rp[sl]")
# Novelty scores
GT@meta.data<-GT@meta.data |>dplyr::rename(nUMI = nCount_RNA, nGene = nFeature_RNA) |>dplyr::mutate(log10GenesPerUMI = log10(nGene) / log10(nUMI))
# Adding one column for sample
GT@meta.data$Sample<- sapply(strsplit(rownames(GT@meta.data), "_"), `[` , 1)
View(GT@meta.data)

#' 
## ------------------------------------------------------------------------------------------------------------------
saveRDS(GT, file = "deconvolution/combined_unfiltered_GT.rds")

#' 
## ----Quality control-----------------------------------------------------------------------------------------------
pdf("deconvolution/QC/QC_control_features.pdf", width = 15, height = 8)
VlnPlot(GT, features = c("nUMI", "nGene", "pMito", "pRibo", "log10GenesPerUMI"), pt.size = 0, group.by = 'Sample', ncol = 4) + theme(legend.position = "none")
invisible(dev.off())

#' 
## ------------------------------------------------------------------------------------------------------------------
#GT<- readRDS("output/combined_unfiltered_GT.rds")
# Filtering
quantile(GT$nUMI,  probs = c(0.001, 0.01, 0.05, 0.5, 0.95, 0.99, 0.999), na.rm = TRUE)
quantile(GT$nGene, probs = c(0.001, 0.01, 0.05, 0.5, 0.95, 0.99, 0.999), na.rm = TRUE)
quantile(GT$pMito, probs = c(0.50, 0.80, 0.90, 0.95, 0.99), na.rm = TRUE)
quantile(GT$log10GenesPerUMI, probs = c(0.01, 0.05, 0.1, 0.5), na.rm = TRUE)

GT_filtered<-subset(GT, nUMI > 300 & nUMI < 10000 & nGene > 250 & pMito <= 13 & log10GenesPerUMI >= 0.93)
saveRDS(GT_filtered, file = "deconvolution/combined_filtered_GT.rds")

#' 
## ------------------------------------------------------------------------------------------------------------------
pdf("deconvolution//QC/QC_after_control_features.pdf", width = 15, height = 8)
VlnPlot(GT_filtered, features = c("nUMI", "nGene", "pMito", "pRibo", "log10GenesPerUMI"), pt.size = 0,group.by = 'Sample', ncol = 4) + theme(legend.position = "none")
invisible(dev.off())

FeatureScatter(GT_filtered, feature1 = "nUMI", feature2 = "nGene")
FeatureScatter(GT_filtered, feature1 = "nUMI", feature2 = "pMito")

#' 
#' #---- Workflow
## ----non integration workflow--------------------------------------------------------------------------------------
#load("MouseCellCycleGenes.rda")
# SCTransform
GT_filtered <- SCTransform(GT_filtered,vars.to.regress = c("pMito", "nUMI"),   verbose = FALSE)
# PCA
GT_filtered <- RunPCA(GT_filtered, verbose = FALSE)
ElbowPlot(GT_filtered, ndims = 50)
VizDimLoadings(GT_filtered, dims = 1:3, reduction = "pca")
DimHeatmap(GT_filtered, dims = 1:15, cells = 500, balanced = TRUE)

GT_filtered <- FindNeighbors(GT_filtered, reduction = "pca",dims = 1:35,nn.eps = 0.5)
GT_filtered <- FindClusters(GT_filtered, resolution = seq(0.1, 1.2, by = 0.1),algorithm = 1,n.iter = 1000)
GT_filtered <- RunUMAP(GT_filtered, dims = 1:35, reduction = "pca", seed.use = 7081998)

DefaultAssay(GT_filtered) <- "RNA"
GT_filtered <- NormalizeData(object = GT_filtered,normalization.method = "LogNormalize", scale.factor = 10000)

DimPlot(GT_filtered, label = T, group.by = "SCT_snn_res.0.5", split.by = "Sample")
FeaturePlot(GT_filtered, label = T, features = "Mki67", cols = c("gray90", "purple4"))

#' 
#' #--- Doublets finding
suppressPackageStartupMessages({
;library(tidyverse);library(ggrepel);library(emmeans)
;library(SingleCellExperiment);library(scater);library(BiocParallel)
;library(ggpubr);library(speckle);library(magrittr)
;library(broom);library(muscat);library(Seurat);library(clustree)
;library(leiden);library(data.table);library(cowplot)
;library(scDblFinder);library(BiocSingular);library(scds)
})

# SingleCellExperiment 
sce <- as.SingleCellExperiment(GT_filtered)

cat("Start running the doubleting detection.....")
set.seed(7081998)
sce <- scDblFinder(sce,samples="Sample", BPPARAM=BiocParallel::SnowParam(workers = 8),nfeatures = 3000,dims = 35, dbr.sd = 1)
cat("Finished..!")
## Adding doublet identification on the metadata
GT_filtered@meta.data$Doublets <- sce$scDblFinder.class

## Subsetting the doublets from singlets 
seuObject_slim_nodoub <- subset(GT_filtered, subset = Doublets == "singlet")

table(sce$scDblFinder.class) ## number of singlets and doublets 
table(sce$scDblFinder.class, sce$Sample) ## number of singlets and doublets per sample
round(100 * prop.table(table(sce$scDblFinder.class, sce$Sample), margin = 2), 2)

#' 
## ----Phase annotation----------------------------------------------------------------------------------------------
load("MouseCellCycleGenes.rda") 
seuObject_slim_nodoub<-CellCycleScoring(seuObject_slim_nodoub,s.features = s_genes,g2m.features = g2m_genes,set.ident = TRUE,nbin = 24)
DimPlot(seuObject_slim_nodoub, group.by = "Phase")

#' 
## ----Testing markers and clustering--------------------------------------------------------------------------------
DimPlot(GT_filtered, label = T, group.by = "SCT_snn_res.0.9", split.by = "Treatment")
DimPlot(seuObject_slim_nodoub, label = T, group.by = "SCT_snn_res.0.9", split.by = "Treatment")

Idents(seuObject_slim_nodoub)<- seuObject_slim_nodoub$SCT_snn_res.0.9
Idents(seuObject_slim_nodoub)<- seuObject_slim_nodoub$SCT_snn_res.0.8
FeaturePlot(seuObject_slim_nodoub, label = T, features = "Top2a", cols = c("gray90", "purple4"),
            min.cutoff = "q10", max.cutoff = "q95")
DimPlot(seuObject_slim_nodoub, label = T, group.by = "SCT_snn_res.0.8")

#' 
## ----Markers-------------------------------------------------------------------------------------------------------
### First running markers for exploration and then for celltype2
## Exploration is for SCT_snn_res.0.8
### Running presto 
## Identity resolution I want 
testing_1 <- seuObject_slim_nodoub
Idents(testing_1) <- "SCT_snn_res.0.8"
DefaultAssay(testing_1)<- "RNA"
## Running Wilcoxau analysis 
all_markers_clustID <- presto::wilcoxauc(testing_1, group_by = 'SCT_snn_res.0.8', assay = 'data')
unique(all_markers_clustID$group)

all_markers_clustID$group <- paste("Cluster",all_markers_clustID$group,sep="_")
all_markers.Sign <- all_markers_clustID %>% dplyr::filter(padj < 0.05, logFC > 0)

top20 <- presto::top_markers(all_markers.Sign,n = 20,auc_min = 0.5, pval_max = 0.05)

openxlsx::write.xlsx(all_markers.Sign,
                     file = "deconvolution/Markers/SCT.0.8_PrestoByCluster_Filteredmarkers_padjLT05_logfcGT0.xlsx",
                     colNames = TRUE,rowNames = FALSE,borders = "columns",sheetName="Markers")
openxlsx::write.xlsx(top20,
                     file = "deconvolution/Markers/SCT.0.8_NOBOUD_Mice_dataset_top20_0.2.xlsx",
                     colNames = TRUE,rowNames = FALSE,borders = "columns",sheetName="Markers")

#' 
## ----General Annotation--------------------------------------------------------------------------------------------
Idents(seuObject_slim_nodoub)<- seuObject_slim_nodoub$SCT_snn_res.0.8
unique(Idents(seuObject_slim_nodoub))
DefaultAssay(seuObject_slim_nodoub)<- "RNA"
cluster_annotations <- c(
  "0"  = "Differentiating CGNP-like tumor cells",
  "1"  = "Cycling tumor cells",
  "2"  = "Early Differentiating CGNP-like cells",
  "3"  = "Granule Cell Precursors",
  "4"  = "Granule Cell Precursors",
  "5"  = "Granule Cell Precursors",
  "6"  = "Cycling tumor cells",
  "7"  = "Cycling tumor cells",
  "8"  = "Early Differentiating CGNP-like cells",
  "9"  = "Stressed/Cycling GCPs",
  "10" = "Oligodendrocyte Precursor Cells",
  "11" = "CGNs and CGN-differentiated tumor cells",
  "12" = "Pvalb+ Tumor Subpopulation",
  "13" = "Astrocytes",
  "14" = "Myeloid cells",
  "15" = "Fibroblasts",
  "16" = "Endothelial Cells"
)
Idents(seuObject_slim_nodoub)<-seuObject_slim_nodoub$SCT_snn_res.0.8
seuObject_slim_nodoub@meta.data$celltype<- cluster_annotations[as.character(Idents(seuObject_slim_nodoub))]

DimPlot(seuObject_slim_nodoub, group.by = "celltype")
DimPlot(seuObject_slim_nodoub, group.by = "Phase")
DimPlot(seuObject_slim_nodoub, split.by = "Treatment")
DimPlot(seuObject_slim_nodoub, group.by = "SCT_snn_res.0.8", label = T)

save(seuObject_slim_nodoub, file = "deconvolution/seuObject_slim_nodoub.RData")

#' 
## ----Vizualization-------------------------------------------------------------------------------------------------
## PCA cladogram
Idents(seuObject_slim_nodoub)<- seuObject_slim_nodoub$SCT_snn_res.0.8
seuObject_slim_nodoub<-BuildClusterTree(seuObject_slim_nodoub, dims = 1:35)  
PlotClusterTree(seuObject_slim_nodoub)

# General umap
Idents(seuObject_slim_nodoub)<- seuObject_slim_nodoub$SCT_snn_res.0.8
nClust <- length(unique(Idents(seuObject_slim_nodoub)))
colCls <- colorRampPalette(brewer.pal(n = 10, name = "Paired"))(nClust)
pdf("deconvolution/General_umap.pdf", width = 8, height = 6)
DimPlot(seuObject_slim_nodoub, reduction = "umap", pt.size = 0.01, label = T, label.size = 5,cols = colCls) +  coord_fixed()
invisible(dev.off())

### by celltype
Idents(seuObject_slim_nodoub)<- seuObject_slim_nodoub$celltype
nClust <- length(unique(Idents(seuObject_slim_nodoub)))
colCls <- colorRampPalette(brewer.pal(n = 10, name = "Paired"))(nClust)
pdf("deconvolution/celltype_General_umapv2.pdf", width = 8, height = 6)
DimPlot(seuObject_slim_nodoub, reduction = "umap", pt.size = 0.01, label = F, label.size = 5,cols = colCls) +  coord_fixed()
invisible(dev.off())

### Proportion by percentage
# Theme
plotTheme <- theme(
  text = element_text(size = 18),
  axis.text.x = element_text(angle = 45, hjust = 1),
  axis.text.y = element_text(size = 16),
  legend.title = element_blank(),
  panel.grid = element_blank(),
  panel.background = element_blank(),
  axis.line = element_line(color = "black")
)

ggData = data.frame(prop.table(table(seuObject_slim_nodoub$SCT_snn_res.0.8, seuObject_slim_nodoub$Treatment), margin = 2))
colnames(ggData) = c("cluster", "sample", "value")
p1 <- ggplot(ggData, aes(sample, value, fill = cluster)) +
  geom_col() + xlab("Sample") + ylab("Proportion of Cells (%)") + scale_fill_manual(values = colCls)+ plotTheme + coord_flip()
p1
ggsave(p1, width = 10, height = 5, filename = "deconvolution/Treatment_proportions.jpg")

### Proportion by sample
ggData = data.frame(prop.table(table(seuObject_slim_nodoub$SCT_snn_res.0.8, seuObject_slim_nodoub$Sample), margin = 2))
colnames(ggData) = c("cluster", "sample", "value")
p1 <- ggplot(ggData, aes(sample, value, fill = cluster)) +
  geom_col() + xlab("Sample") + ylab("Proportion of Cells (%)") + scale_fill_manual(values = colCls)+ plotTheme + coord_flip()
p1
ggsave(p1,  width = 10, height = 5, filename = "deconvolution/Sample_proportions.jpg")



#' 
#' 
#' #--------------------------------------
#' ## Section 2: Analysis on TME
## ----TME workflow--------------------------------------------------------------------------------------------------
DefaultAssay(seuObject_slim_nodoub)<- "RNA"
Idents(seuObject_slim_nodoub)<- seuObject_slim_nodoub$celltype
unique(seuObject_slim_nodoub$celltype)

# Subset TME
TME<- subset(seuObject_slim_nodoub, idents = c("Astrocytes", "Oligodendrocyte Precursor Cells", "Endothelial Cells",
                                               "Fibroblasts", "Myeloid cells"))
DimPlot(TME, split.by = "Treatment")

### source and re do analysis 
source("output/utils/Utils.R")
TME<- processing_seurat_sctransform(TME, vars_to_regress = c("nUMI", "pMito"), 
                                    npcs = 35, res = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1, 1.1, 1.2))

DefaultAssay(TME) <- "RNA"
TME <- NormalizeData(object = TME,normalization.method = "LogNormalize", scale.factor = 10000)

DimPlot(TME, group.by = "SCT_snn_res.0.6", label = T)
FeaturePlot(TME, features = "nUMI", split.by = "Treatment", label = T)
VlnPlot(TME, features = "nUMI")

### First running markers for exploration.
## Exploration is for SCT_snn_res.0.6
## Markers
testing_1<- TME
Idents(testing_1) <- "SCT_snn_res.0.6"
## Running Wilcoxau analysis 
all_markers_clustID <- presto::wilcoxauc(testing_1, group_by = 'SCT_snn_res.0.6', assay = 'data')
unique(all_markers_clustID$group)

all_markers_clustID$group <- paste("Cluster",all_markers_clustID$group,sep="_")
all_markers.Sign <- all_markers_clustID %>%dplyr::filter(padj < 0.05, logFC > 0)
top20 <- presto::top_markers(all_markers.Sign,n = 20,auc_min = 0.5, pval_max = 0.05)

openxlsx::write.xlsx(all_markers.Sign,
                     file = "deconvolution/TME/SCT_snn_res.0.6_PrestoByCluster_Filteredmarkers_padjLT05_logfcGT0.xlsx",
                     colNames = TRUE,rowNames = FALSE,borders = "columns",sheetName="Markers")
openxlsx::write.xlsx(top20,
                     file = "deconvolution/TME/SCT_snn_res.0.6_NOBOUD_Mice_dataset_top20_0.2.xlsx",
                     colNames = TRUE,rowNames = FALSE,borders = "columns",sheetName="Markers")

DimPlot(TME, group.by = "SCT_snn_res.0.6", label = T)

#' 
## ----TME annotation and viz----------------------------------------------------------------------------------------
cluster_annotations <- c(
  "0"  = "OPCs-neuronal",
  "1"  = "Myeloid cells",
  "2"  = "Astrocytes",
  "3"  = "Fibroblasts CAFs",
  "4"  = "Endothelial cells",
  "5"  = "Mature oligodendrocytes"
)
Idents(TME)<-TME$SCT_snn_res.0.6
TME@meta.data$celltype2<- cluster_annotations[as.character(Idents(TME))]
DimPlot(TME, group.by = "celltype2", label = T)
save(TME, file = "deconvolution/TME/TME.RData")

## Viz 
### General umap
### Visualization object 
Idents(TME)<- TME$celltype2
nClust <- length(unique(Idents(TME)))# Setup color palette
colCls <- colorRampPalette(brewer.pal(n = 10, name = "Paired"))(nClust)
pdf("deconvolution/TME/General_umap.pdf", width = 8, height = 6)
DimPlot(TME, reduction = "umap", pt.size = 0.8, label = F, label.size = 4,cols = colCls) +  coord_fixed()
invisible(dev.off())

# Theme
plotTheme <- theme(
  text = element_text(size = 18),
  axis.text.x = element_text(angle = 45, hjust = 1),
  axis.text.y = element_text(size = 16),
  legend.title = element_blank(),
  panel.grid = element_blank(),
  panel.background = element_blank(),
  axis.line = element_line(color = "black")
)

ggData = data.frame(prop.table(table(TME$celltype2, TME$Treatment), margin = 2))
colnames(ggData) = c("cluster", "sample", "value")
p1 <- ggplot(ggData, aes(sample, value, fill = cluster)) + geom_col() + xlab("Sample") + ylab("Proportion of Cells (%)") +
  scale_fill_manual(values = colCls)+ plotTheme + coord_flip()
ggsave(p1, width = 10, height = 5, filename = "deconvolution/TME/Treatment_proportions.jpg")

# Proportion by sample
ggData = data.frame(prop.table(table(TME$celltype2, TME$Sample), margin = 2))
colnames(ggData) = c("cluster", "sample", "value")
p1 <- ggplot(ggData, aes(sample, value, fill = cluster)) +geom_col() + xlab("Sample") + ylab("Proportion of Cells (%)") +
  scale_fill_manual(values = colCls)+ plotTheme + coord_flip()
ggsave(p1, width = 10, height = 5, filename = "deconvolution/TME/Sample_proportions.jpg")


#' 
# Session Information

# R version 4.5.1 (2025-06-13 ucrt)
# Platform: x86_64-w64-mingw32/x64
# Running under: Windows 11 x64 (build 26200)
# 
# Matrix products: default
#   LAPACK version 3.12.1
# 
# locale:
# [1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8   
# [3] LC_MONETARY=English_United States.utf8 LC_NUMERIC=C                          
# [5] LC_TIME=English_United States.utf8    
# 
# time zone: America/New_York
# tzcode source: internal
# 
# attached base packages:
# [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#  [1] scds_1.24.0                 BiocSingular_1.26.1         scDblFinder_1.22.0          leiden_0.4.3.1             
#  [5] clustree_0.5.1              ggraph_2.2.2                muscat_1.22.0               broom_1.0.12               
#  [9] magrittr_2.0.4              speckle_1.8.0               ggpubr_0.6.3                BiocParallel_1.44.0        
# [13] scater_1.36.0               emmeans_2.0.3               ggrepel_0.9.8               sessioninfo_1.2.3          
# [17] RColorBrewer_1.1-3          presto_1.0.0                data.table_1.18.2.1         Rcpp_1.1.1                 
# [21] scuttle_1.18.0              SeuratWrappers_0.4.0        knitr_1.51                  openxlsx_4.2.8.1           
# [25] pheatmap_1.0.13             RCurl_1.98-1.18             cowplot_1.2.0               scales_1.4.0               
# [29] Matrix_1.7-3                lubridate_1.9.5             forcats_1.0.1               stringr_1.6.0              
# [33] dplyr_1.2.1                 purrr_1.2.2                 readr_2.2.0                 tidyr_1.3.2                
# [37] tibble_3.3.1                ggplot2_4.0.2               tidyverse_2.0.0             Seurat_5.4.0               
# [41] SeuratObject_5.3.0          sp_2.2-1                    SingleCellExperiment_1.30.1 SummarizedExperiment_1.38.1
# [45] Biobase_2.68.0              GenomicRanges_1.60.0        GenomeInfoDb_1.44.3         IRanges_2.44.0             
# [49] S4Vectors_0.48.0            BiocGenerics_0.54.1         generics_0.1.4              MatrixGenerics_1.20.0      
# [53] matrixStats_1.5.0          
# 
# loaded via a namespace (and not attached):
#   [1] R.methodsS3_1.8.2        dichromat_2.0-0.1        progress_1.2.3           goftest_1.2-3           
#   [5] Biostrings_2.78.0        vctrs_0.7.2              spatstat.random_3.4-5    digest_0.6.39           
#   [9] png_0.1-9                corpcor_1.6.10           shape_1.4.6.1            deldir_2.0-4            
#  [13] parallelly_1.46.1        MASS_7.3-65              reshape2_1.4.5           httpuv_1.6.17           
#  [17] foreach_1.5.2            withr_3.0.2              xfun_0.55                survival_3.8-3          
#  [21] memoise_2.0.1            ggbeeswarm_0.7.3         Seqinfo_1.0.0            zoo_1.8-15              
#  [25] GlobalOptions_0.1.4      gtools_3.9.5             pbapply_1.7-4            R.oo_1.27.1             
#  [29] Formula_1.2-5            prettyunits_1.2.0        promises_1.5.0           otel_0.2.0              
#  [33] httr_1.4.8               restfulr_0.0.16          rstatix_0.7.3            globals_0.19.1          
#  [37] fitdistrplus_1.2-6       rstudioapi_0.18.0        UCSC.utils_1.4.0         miniUI_0.1.2            
#  [41] curl_7.0.0               ScaledMatrix_1.18.0      polyclip_1.10-7          GenomeInfoDbData_1.2.14 
#  [45] SparseArray_1.10.8       xtable_1.8-8             doParallel_1.0.17        evaluate_1.0.5          
#  [49] S4Arrays_1.10.1          hms_1.1.4                irlba_2.3.7              colorspace_2.1-2        
#  [53] ROCR_1.0-12              reticulate_1.46.0        spatstat.data_3.1-9      lmtest_0.9-40           
#  [57] later_1.4.8              viridis_0.6.5            lattice_0.22-7           spatstat.geom_3.7-3     
#  [61] future.apply_1.20.2      XML_3.99-0.23            scattermore_1.2          RcppAnnoy_0.0.23        
#  [65] pillar_1.11.1            nlme_3.1-168             iterators_1.0.14         caTools_1.18.3          
#  [69] compiler_4.5.1           beachmat_2.26.0          RSpectra_0.16-2          stringi_1.8.7           
#  [73] tensor_1.5.1             minqa_1.2.8              GenomicAlignments_1.44.0 plyr_1.8.9              
#  [77] BiocIO_1.18.0            crayon_1.5.3             abind_1.4-8              blme_1.0-7              
#  [81] locfit_1.5-9.12          graphlayouts_1.2.3       sandwich_3.1-1           codetools_0.2-20        
#  [85] GetoptLong_1.1.1         plotly_4.12.0            remaCor_0.0.20           mime_0.13               
#  [89] splines_4.5.1            circlize_0.4.18          fastDummies_1.7.5        here_1.0.2              
#  [93] clue_0.3-68              lme4_2.0-1               listenv_0.10.1           Rdpack_2.6.6            
#  [97] ggsignif_0.6.4           estimability_1.5.1       statmod_1.5.1            tzdb_0.5.0              
# [101] fANCOVA_0.6-1            tweenr_2.0.3             pkgconfig_2.0.3          tools_4.5.1             
# [105] cachem_1.1.0             RhpcBLASctl_0.23-42      rbibutils_2.4.1          viridisLite_0.4.3       
# [109] numDeriv_2016.8-1.1      fastmap_1.2.0            grid_4.5.1               ica_1.0-3               
# [113] Rsamtools_2.24.1         patchwork_1.3.2          coda_0.19-4.1            BiocManager_1.30.27     
# [117] dotCall64_1.2            carData_3.0-6            RANN_2.6.2               farver_2.1.2            
# [121] reformulas_0.4.4         aod_1.3.3                tidygraph_1.3.1          mgcv_1.9-3              
# [125] yaml_2.3.12              rtracklayer_1.68.0       cli_3.6.6                lifecycle_1.0.5         
# [129] uwot_0.2.4               glmmTMB_1.1.14           mvtnorm_1.3-6            bluster_1.18.0          
# [133] backports_1.5.1          timechange_0.4.0         gtable_0.3.6             rjson_0.2.23            
# [137] ggridges_0.5.7           progressr_0.19.0         pROC_1.19.0.1            parallel_4.5.1          
# [141] limma_3.66.0             jsonlite_2.0.0           edgeR_4.8.2              RcppHNSW_0.6.0          
# [145] bitops_1.0-9             xgboost_3.2.1.1          Rtsne_0.17               spatstat.utils_3.2-2    
# [149] BiocNeighbors_2.2.0      zip_2.3.3                metapod_1.16.0           dqrng_0.4.1             
# [153] spatstat.univar_3.1-7    R.utils_2.13.0           pbkrtest_0.5.5           lazyeval_0.2.3          
# [157] shiny_1.13.0             htmltools_0.5.9          sctransform_0.4.3        rappdirs_0.3.4          
# [161] glue_1.8.0               spam_2.11-3              XVector_0.50.0           rprojroot_2.1.1         
# [165] scran_1.36.0             gridExtra_2.3            EnvStats_3.1.0           boot_1.3-31             
# [169] igraph_2.2.3             variancePartition_1.38.1 TMB_1.9.21               R6_2.6.1                
# [173] DESeq2_1.48.2            gplots_3.3.0             cluster_2.1.8.2          nloptr_2.2.1            
# [177] DelayedArray_0.36.0      tidyselect_1.2.1         vipor_0.4.7              ggforce_0.5.0           
# [181] car_3.1-5                future_1.70.0            rsvd_1.0.5               KernSmooth_2.23-26      
# [185] S7_0.2.1                 htmlwidgets_1.6.4        ComplexHeatmap_2.26.0    rlang_1.1.7             
# [189] spatstat.sparse_3.1-0    spatstat.explore_3.8-0   lmerTest_3.2-1           remotes_2.5.0           
# [193] beeswarm_0.4.0          
