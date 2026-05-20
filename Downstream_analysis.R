#' ---
#' title: "Downstream analysis on DEGs from DESEQ model"
#' output: html_document
#' date: "2025-11-25"
#' ---
#' 
## ----setup, include=FALSE------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

#' 
## ----Libraries-----------------------------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(ggplot2);library(ggrepel);library(ggpubr);
  library(data.table);library(RColorBrewer);library(tidyverse);
  library(preprocessCore);library(future.apply);library(DESeq2);
  library(pheatmap);library(sva);library(viridis);
  library(limma);library(emmeans);library(broom);
  library(janitor);library(tidyplots);library(dplyr);
  library(writexl);library(scToppR);library(fgsea);
  library(clusterProfiler);library(readxl);
  library(ComplexHeatmap);library(sessioninfo);
  library(org.Mm.eg.db)})

#' 
## ----Loading Gene dataframes and preprocess------------------------------------------------------------------------

Combo <- read_xlsx("edgeR_new_exploration/deg/Combo_DMSO/CTX_FullTab_shrink.xlsx") # Combo vs DMSO
BMS <- read_xlsx("edgeR_new_exploration/deg/BMS_DMSO/CTX_FullTab_shrink.xlsx") # BMS vs DMSO
Ribo <- read_xlsx("edgeR_new_exploration/deg/Ribo_DMSO/CTX_FullTab_shrink.xlsx") # Ribo vs DMSO

Combo_BMS<- read_xlsx("edgeR_new_exploration/deg/Combo_BMS/CTX_FullTab_shrink.xlsx") # Combo vs BMS
Ribo_BMS<- read_xlsx("edgeR_new_exploration/deg/Ribo_BMS/CTX_FullTab_shrink.xlsx") # Ribo vs BMS

Combo_Ribo<- read_xlsx("edgeR_new_exploration/deg/Combo_Ribo/CTX_FullTab_shrink.xlsx") # Combo vs Ribo

Unique_Combo<-read_xlsx("edgeR_new_exploration/deg/Unique_combo/Combo_unique_CTX_DGE_shrink.xlsx") # Unique Combo genes
Unique_BMS<- read_xlsx("edgeR_new_exploration/deg/Unique_BMS/unique_BMS_DEGs.xlsx") # Unique BMS genes
Unique_Ribo<- read_xlsx("edgeR_new_exploration/deg/Unique_Ribo/unique_Ribo_DEGs.xlsx") # Unique Ribo genes

## Creating tmps
Combo <- Combo %>%
  mutate(LOG = -log10(padj), ABS = abs(log2FoldChange)) %>% mutate(Threshold = if_else(padj < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>% mutate(Direction = case_when(log2FoldChange > 0.3 & padj < 0.05 ~ "UpReg", log2FoldChange < -0.3 & padj < 0.05 ~ "DownReg"))

BMS <- BMS %>%
  mutate(LOG = -log10(padj), ABS = abs(log2FoldChange)) %>% mutate(Threshold = if_else(padj < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>% mutate(Direction = case_when(log2FoldChange > 0.3 & padj < 0.05 ~ "UpReg", log2FoldChange < -0.3 & padj < 0.05 ~ "DownReg"))

Ribo <- Ribo %>%
  mutate(LOG = -log10(padj), ABS = abs(log2FoldChange)) %>% mutate(Threshold = if_else(padj < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>% mutate(Direction = case_when(log2FoldChange > 0.3 & padj < 0.05 ~ "UpReg", log2FoldChange < -0.3 & padj < 0.05 ~ "DownReg"))

Combo_BMS <- Combo_BMS %>%
  mutate(LOG = -log10(padj), ABS = abs(log2FoldChange)) %>% mutate(Threshold = if_else(padj < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>% mutate(Direction = case_when(log2FoldChange > 0.3 & padj < 0.05 ~ "UpReg", log2FoldChange < -0.3 & padj < 0.05 ~ "DownReg"))

Ribo_BMS <- Ribo_BMS %>%
  mutate(LOG = -log10(padj), ABS = abs(log2FoldChange)) %>% mutate(Threshold = if_else(padj < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>% mutate(Direction = case_when(log2FoldChange > 0.3 & padj < 0.05 ~ "UpReg", log2FoldChange < -0.3 & padj < 0.05 ~ "DownReg"))

Combo_Ribo <- Combo_Ribo %>%
  mutate(LOG = -log10(padj), ABS = abs(log2FoldChange)) %>% mutate(Threshold = if_else(padj < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>% mutate(Direction = case_when(log2FoldChange > 0.3 & padj < 0.05 ~ "UpReg", log2FoldChange < -0.3 & padj < 0.05 ~ "DownReg"))

Unique_Combo <- Unique_Combo %>%
  mutate(LOG = -log10(padj), ABS = abs(log2FoldChange)) %>% mutate(Threshold = if_else(padj < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>% mutate(Direction = case_when(log2FoldChange > 0.3 & padj < 0.05 ~ "UpReg", log2FoldChange < -0.3 & padj < 0.05 ~ "DownReg"))

Unique_BMS <- Unique_BMS %>%
  mutate(LOG = -log10(padj), ABS = abs(log2FoldChange)) %>% mutate(Threshold = if_else(padj < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>% mutate(Direction = case_when(log2FoldChange > 0.3 & padj < 0.05 ~ "UpReg", log2FoldChange < -0.3 & padj < 0.05 ~ "DownReg"))

Unique_Ribo <- Unique_Ribo %>%
  mutate(LOG = -log10(padj), ABS = abs(log2FoldChange)) %>% mutate(Threshold = if_else(padj < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>% mutate(Direction = case_when(log2FoldChange > 0.3 & padj < 0.05 ~ "UpReg", log2FoldChange < -0.3 & padj < 0.05 ~ "DownReg"))

#' 
#' #--------------------------------------
#' # Gene Ontology on each comparision
#--------------------------------------
# Gene Ontology on Combo vs DMSO DEGs
tmp<- Combo |> filter(Direction %in% c("UpReg", "DownReg"))
toppData <-  toppFun(
       tmp,topp_categories = NULL,cluster_col = "Direction",gene_col = "Gene", p_val_col = "padj",
       logFC_col = "log2FoldChange", pval_cutoff = 0.05, min_genes = 10, max_genes = 500, max_results = 50)
toppPlot(
       toppData,category = "GeneOntologyMolecularFunction",num_terms = 10,
       p_val_adj = "BH",p_val_display = "log",save = TRUE,
       save_dir = "edgeR_new_exploration/deg/Combo_DMSO/GO/pseudobulk_moderate_GO",
       width = 5, height = 6)
toppPlot(
       toppData,category = "GeneOntologyBiologicalProcess",num_terms = 10,
       p_val_adj = "BH",p_val_display = "log",save = TRUE,
       save_dir = "edgeR_new_exploration/deg/Combo_DMSO/GO/pseudobulk_moderate_GO",
       width = 5, height = 6)

GO<- toppData |> as.data.frame()
write_xlsx(GO, "edgeR_new_exploration/deg/Combo_DMSO/GO/pseudobulk_moderate_GO/GO.xlsx")
save(toppData, file = "edgeR_new_exploration/deg/Combo_DMSO/GO/pseudobulk_moderate_GO/TopData.RData")

#--------------------------------------
# Gene Ontology on BMS vs DMSO DEGs
tmp<- BMS |> filter(Direction %in% c("UpReg", "DownReg"))
toppData <-  toppFun(
       tmp,topp_categories = NULL,cluster_col = "Direction",gene_col = "Gene", p_val_col = "padj",
       logFC_col = "log2FoldChange", pval_cutoff = 0.05, min_genes = 10, max_genes = 500, max_results = 50)
toppPlot(
       toppData,category = "GeneOntologyMolecularFunction",num_terms = 10,
       p_val_adj = "BH",p_val_display = "log",save = TRUE,
       save_dir = "edgeR_new_exploration/deg/BMS_DMSO/GO/pseudobulk_moderate_GO",
       width = 5, height = 6)
toppPlot(
       toppData,category = "GeneOntologyBiologicalProcess",num_terms = 10,
       p_val_adj = "BH",p_val_display = "log",save = TRUE,
       save_dir = "edgeR_new_exploration/deg/BMS_DMSO/GO/pseudobulk_moderate_GO",
       width = 5, height = 6)

GO<- toppData |> as.data.frame()
write_xlsx(GO, "edgeR_new_exploration/deg/BMS_DMSO/GO/pseudobulk_moderate_GO/GO.xlsx")
save(toppData, file = "edgeR_new_exploration/deg/BMS_DMSO/GO/pseudobulk_moderate_GO/TopData.RData")

#--------------------------------------
# Gene Ontology on Ribo vs DMSO DEGs
tmp<- Ribo |> filter(Direction %in% c("UpReg", "DownReg"))
toppData <-  toppFun(
       tmp,topp_categories = NULL,cluster_col = "Direction",gene_col = "Gene", p_val_col = "padj",
       logFC_col = "log2FoldChange", pval_cutoff = 0.05, min_genes = 10, max_genes = 500, max_results = 50)
toppPlot(
       toppData,category = "GeneOntologyMolecularFunction",num_terms = 10,
       p_val_adj = "BH",p_val_display = "log",save = TRUE,
       save_dir = "edgeR_new_exploration/deg/Ribo_DMSO/GO/pseudobulk_moderate_GO",
       width = 5, height = 6)
toppPlot(
       toppData,category = "GeneOntologyBiologicalProcess",num_terms = 10,
       p_val_adj = "BH",p_val_display = "log",save = TRUE,
       save_dir = "edgeR_new_exploration/deg/Ribo_DMSO/GO/pseudobulk_moderate_GO",
       width = 5, height = 6)

GO<- toppData |> as.data.frame()
write_xlsx(GO, "edgeR_new_exploration/deg/Ribo_DMSO/GO/pseudobulk_moderate_GO/GO.xlsx")
save(toppData, file = "edgeR_new_exploration/deg/Ribo_DMSO/GO/pseudobulk_moderate_GO/TopData.RData")

#--------------------------------------
# Gene Ontology on Combo vs BMS DEGs
tmp<- Combo_BMS |> filter(Direction %in% c("UpReg", "DownReg"))
toppData <-  toppFun(
       tmp,topp_categories = NULL,cluster_col = "Direction",gene_col = "Gene", p_val_col = "padj",
       logFC_col = "log2FoldChange", pval_cutoff = 0.05, min_genes = 10, max_genes = 500, max_results = 50)
toppPlot(
       toppData,category = "GeneOntologyMolecularFunction",num_terms = 10,
       p_val_adj = "BH",p_val_display = "log",save = TRUE,
       save_dir = "edgeR_new_exploration/deg/Combo_BMS/GO/pseudobulk_moderate_GO",
       width = 5, height = 6)
toppPlot(
       toppData,category = "GeneOntologyBiologicalProcess",num_terms = 10,
       p_val_adj = "BH",p_val_display = "log",save = TRUE,
       save_dir = "edgeR_new_exploration/deg/Combo_BMS/GO/pseudobulk_moderate_GO",
       width = 5, height = 6)

GO<- toppData |> as.data.frame()
write_xlsx(GO, "edgeR_new_exploration/deg/Combo_BMS/GO/pseudobulk_moderate_GO/GO.xlsx")
save(toppData, file = "edgeR_new_exploration/deg/Combo_BMS/GO/pseudobulk_moderate_GO/TopData.RData")

#--------------------------------------
# Gene Ontology on Ribo vs BMS DEGs
tmp<- Ribo_BMS |> filter(Direction %in% c("UpReg", "DownReg"))
toppData <-  toppFun(
       tmp,topp_categories = NULL,cluster_col = "Direction",gene_col = "Gene", p_val_col = "padj",
       logFC_col = "log2FoldChange", pval_cutoff = 0.05, min_genes = 10, max_genes = 500, max_results = 50)
toppPlot(
       toppData,category = "GeneOntologyMolecularFunction",num_terms = 10,
       p_val_adj = "BH",p_val_display = "log",save = TRUE,
       save_dir = "edgeR_new_exploration/deg/Ribo_BMS/GO/pseudobulk_moderate_GO",
       width = 5, height = 6)
toppPlot(
       toppData,category = "GeneOntologyBiologicalProcess",num_terms = 10,
       p_val_adj = "BH",p_val_display = "log",save = TRUE,
       save_dir = "edgeR_new_exploration/deg/Ribo_BMS/GO/pseudobulk_moderate_GO",
       width = 5, height = 6)

GO<- toppData |> as.data.frame()
write_xlsx(GO, "edgeR_new_exploration/deg/Ribo_BMS/GO/pseudobulk_moderate_GO/GO.xlsx")
save(toppData, file = "edgeR_new_exploration/deg/Ribo_BMS/GO/pseudobulk_moderate_GO/TopData.RData")

#--------------------------------------
# Gene Ontology on Combo vs Ribo DEGs
tmp<- Combo_Ribo |> filter(Direction %in% c("UpReg", "DownReg"))
toppData <-  toppFun(
       tmp,topp_categories = NULL,cluster_col = "Direction",gene_col = "Gene", p_val_col = "padj",
       logFC_col = "log2FoldChange", pval_cutoff = 0.05, min_genes = 10, max_genes = 500, max_results = 50)
toppPlot(
       toppData,category = "GeneOntologyMolecularFunction",num_terms = 10,
       p_val_adj = "BH",p_val_display = "log",save = TRUE,
       save_dir = "edgeR_new_exploration/deg/Combo_Ribo/GO/pseudobulk_moderate_GO",
       width = 5, height = 6)
toppPlot(
       toppData,category = "GeneOntologyBiologicalProcess",num_terms = 10,
       p_val_adj = "BH",p_val_display = "log",save = TRUE,
       save_dir = "edgeR_new_exploration/deg/Combo_Ribo/GO/pseudobulk_moderate_GO",
       width = 5, height = 6)

GO<- toppData |> as.data.frame()
write_xlsx(GO, "edgeR_new_exploration/deg/Combo_Ribo/GO/pseudobulk_moderate_GO/GO.xlsx")
save(toppData, file = "edgeR_new_exploration/deg/Combo_Ribo/GO/pseudobulk_moderate_GO/TopData.RData")

#--------------------------------------
# Gene Ontology on Unique DEGs from Combo
tmp<- Unique_Combo |> filter(Direction %in% c("UpReg", "DownReg"))
toppData <-  toppFun(
       tmp,topp_categories = NULL,cluster_col = "Direction",gene_col = "Gene", p_val_col = "padj",
       logFC_col = "log2FoldChange", pval_cutoff = 0.05, min_genes = 10, max_genes = 500, max_results = 50)
toppPlot(
       toppData,category = "GeneOntologyMolecularFunction",num_terms = 10,
       p_val_adj = "BH",p_val_display = "log",save = TRUE,
       save_dir = "edgeR_new_exploration/deg/Unique_combo/GO/pseudobulk_moderate_GO",
       width = 5, height = 6)
toppPlot(
       toppData,category = "GeneOntologyBiologicalProcess",num_terms = 10,
       p_val_adj = "BH",p_val_display = "log",save = TRUE,
       save_dir = "edgeR_new_exploration/deg/Unique_combo/GO/pseudobulk_moderate_GO",
       width = 5, height = 6)

GO<- toppData |> as.data.frame()
write_xlsx(GO, "edgeR_new_exploration/deg/Unique_combo/GO/pseudobulk_moderate_GO/GO.xlsx")
save(toppData, file = "edgeR_new_exploration/deg/Unique_combo/GO/pseudobulk_moderate_GO/TopData.RData")

#--------------------------------------
# Gene Ontology on Unique DEGs from BMS
tmp<- Unique_BMS |> filter(Direction %in% c("UpReg", "DownReg"))
toppData <-  toppFun(
       tmp,topp_categories = NULL,cluster_col = "Direction",gene_col = "Gene", p_val_col = "padj",
       logFC_col = "log2FoldChange", pval_cutoff = 0.05, min_genes = 10, max_genes = 500, max_results = 50)
toppPlot(
       toppData,category = "GeneOntologyMolecularFunction",num_terms = 10,
       p_val_adj = "BH",p_val_display = "log",save = TRUE,
       save_dir = "edgeR_new_exploration/deg/Unique_BMS/GO/pseudobulk_moderate_GO",
       width = 5, height = 6)
toppPlot(
       toppData,category = "GeneOntologyBiologicalProcess",num_terms = 10,
       p_val_adj = "BH",p_val_display = "log",save = TRUE,
       save_dir = "edgeR_new_exploration/deg/Unique_BMS/GO/pseudobulk_moderate_GO",
       width = 5, height = 6)

GO<- toppData |> as.data.frame()
write_xlsx(GO, "edgeR_new_exploration/deg/Unique_BMS/GO/pseudobulk_moderate_GO/GO.xlsx")
save(toppData, file = "edgeR_new_exploration/deg/Unique_BMS/GO/pseudobulk_moderate_GO/TopData.RData")

#--------------------------------------
# Gene Ontology on Unique DEGs from Ribo
tmp<- Unique_Ribo |> filter(Direction %in% c("UpReg", "DownReg"))
toppData <-  toppFun(
       tmp,topp_categories = NULL,cluster_col = "Direction",gene_col = "Gene", p_val_col = "padj",
       logFC_col = "log2FoldChange", pval_cutoff = 0.05, min_genes = 10, max_genes = 500, max_results = 50)
toppPlot(
       toppData,category = "GeneOntologyMolecularFunction",num_terms = 10,
       p_val_adj = "BH",p_val_display = "log",save = TRUE,
       save_dir = "edgeR_new_exploration/deg/Unique_Ribo/GO/pseudobulk_moderate_GO",
       width = 5, height = 6)
toppPlot(
       toppData,category = "GeneOntologyBiologicalProcess",num_terms = 10,
       p_val_adj = "BH",p_val_display = "log",save = TRUE,
       save_dir = "edgeR_new_exploration/deg/Unique_Ribo/GO/pseudobulk_moderate_GO",
       width = 5, height = 6)

GO<- toppData |> as.data.frame()
write_xlsx(GO, "edgeR_new_exploration/deg/Unique_Ribo/GO/pseudobulk_moderate_GO/GO.xlsx")
save(toppData, file = "edgeR_new_exploration/deg/Unique_Ribo/GO/pseudobulk_moderate_GO/TopData.RData")

#' 
## ----Gene ontology on AMPK-related pathways from Unique Combo vs DMSO genes----------------------------------------
load("edgeR_new_exploration/deg/Unique_combo/GO/pseudobulk_moderate_GO/TopData.RData")
x <- toppData
x <- x |> filter(Name %in% c("AMP-activated protein kinase activity",
                            "3-phosphoinositide-dependent protein kinase activity",
                            "ribosomal protein S6 kinase activity",
                            "eukaryotic translation initiation factor 2alpha kinase activity",
                            "insulin receptor binding",
                            "ubiquitin-modified protein reader activity"))

toppPlot(x,category = "GeneOntologyMolecularFunction",num_terms = 8,p_val_adj = "BH",
       p_val_display = "log",save = TRUE,
       save_dir = "edgeR_new_exploration/deg/Unique_combo/AMPK_signaling/", width = 6.5,height = 4)

#' 
#' #--------------------------------------
#' # Gene Set Enrichment analysis on Combo, Ribo & BMS vs DMSO
## ------------------------------------------------------------------------------------------------------------------
# Function to run GSEA
run_gsea_per_cluster <- function(df, pathway_list, gene_col= "feature", stat_col = "logFC",
                                 cluster_col = "group", output_dir  = "deg/GSEA/", nperm = 10000, minSize = 10,
                                 maxSize= 500,padj_cutoff = 0.05,label_prefix = "contrast") {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  clusts <- sort(unique(df[[cluster_col]]))
  message("Running GSEA for ", length(clusts), " clusters: ", paste(clusts, collapse = ", "))
  per_cluster <- map(clusts, function(cl) {
    sub <- df %>% filter(.data[[cluster_col]] == cl) %>%
      filter(!is.na(.data[[stat_col]])) %>%
      # If a gene appears multiple times within a cluster, keep the largest |stat|
      arrange(desc(abs(.data[[stat_col]]))) %>% distinct(.data[[gene_col]], .keep_all = TRUE) %>%
      arrange(desc(.data[[stat_col]]))
    stat_vec <- sub[[stat_col]]
    names(stat_vec) <- sub[[gene_col]]
    if (length(stat_vec) < minSize) {
      warning("Cluster ", cl, ": fewer than minSize genes (", length(stat_vec), "). Skipping.")
      return(tibble())
    }
    set.seed(0708)
    fg <- fgsea(pathways = pathway_list, stats = stat_vec, nperm = nperm, minSize = minSize, maxSize= maxSize)
    res <- as_tibble(fg) %>% arrange(padj, desc(NES)) %>% mutate(cluster = cl, label = label_prefix)
    write_csv(res, file.path(output_dir, sprintf("%s_cluster-%s_GSEA.csv", label_prefix, cl)))
    # significant up/down
    up   <- res %>% filter(padj < padj_cutoff, ES > 0)
    down <- res %>% filter(padj < padj_cutoff, ES < 0)
    write_csv(up,   file.path(output_dir, sprintf("%s_cluster-%s_GSEA_up.csv",   label_prefix, cl)))
    write_csv(down, file.path(output_dir, sprintf("%s_cluster-%s_GSEA_down.csv", label_prefix, cl)))
    res
  })
  combined <- bind_rows(per_cluster)
  write_csv(combined, file.path(output_dir, sprintf("%s_ALL_clusters_GSEA.csv", label_prefix)))
  invisible(list(per_cluster = per_cluster, combined = combined))
}

# Loading dataframes and hallmarks
Combo <- read_xlsx("edgeR_new_exploration/deg/Combo_DMSO/CTX_FullTab_shrink.xlsx") # Combo vs DMSO
BMS <- read_xlsx("edgeR_new_exploration/deg/BMS_DMSO/CTX_FullTab_shrink.xlsx") # BMS vs DMSO
Ribo <- read_xlsx("edgeR_new_exploration/deg/Ribo_DMSO/CTX_FullTab_shrink.xlsx") # Ribo vs DMSO
hallmark_pathways <- gmtPathways("edgeR_new_exploration/deg/BMS_DMSO/GSEA/mh.all.v2024.1.Mm.symbols.gmt") # mouse hallmarks

## -- GSEA on Combo vs DMSO
Combo$cluster<- "CombovsDMSO"
res <- run_gsea_per_cluster(df = Combo, pathway_list = hallmark_pathways, gene_col = "Gene",stat_col = "log2FoldChange",
  cluster_col  = "cluster", output_dir   = "edgeR_new_exploration/deg/Combo_DMSO/GSEA/",nperm = 10000,
  label_prefix = "Combo_vs_DMSO"
)
res$combined %>% arrange(padj) %>% head(50)
gsea_output1 <- res$combined |> as.data.frame()
gsea_output1$leadingEdge <- sapply(gsea_output1$leadingEdge, function(x) paste(x, collapse = ", "))
write.csv(gsea_output1, "edgeR_new_exploration/deg/Combo_DMSO/GSEA/full_genes_list_GSEA.csv")
save(res, file = 'edgeR_new_exploration/deg/Combo_DMSO/GSEA/GSEA.RData')

## -- GSEA on BMS vs DMSO
BMS$cluster<- "BMSvsDMSO"
res <- run_gsea_per_cluster(df = BMS, pathway_list = hallmark_pathways, gene_col = "Gene",stat_col = "log2FoldChange",
  cluster_col  = "cluster", output_dir   = "edgeR_new_exploration/deg/BMS_DMSO/GSEA/",nperm = 10000,
  label_prefix = "BMS_vs_DMSO"
)
res$combined %>% arrange(padj) %>% head(50)
gsea_output1 <- res$combined |> as.data.frame()
gsea_output1$leadingEdge <- sapply(gsea_output1$leadingEdge, function(x) paste(x, collapse = ", "))
write.csv(gsea_output1, "edgeR_new_exploration/deg/BMS_DMSO/GSEA/full_genes_list_GSEA.csv")
save(res, file = 'edgeR_new_exploration/deg/BMS_DMSO/GSEA/GSEA.RData')

## -- GSEA on Ribo vs DMSO
Ribo$cluster<- "RibovsDMSO"
res <- run_gsea_per_cluster(df = Ribo, pathway_list = hallmark_pathways, gene_col = "Gene",stat_col = "log2FoldChange",
  cluster_col  = "cluster", output_dir   = "edgeR_new_exploration/deg/Ribo_DMSO/GSEA/",nperm = 10000,
  label_prefix = "Ribo_vs_DMSO"
)
res$combined %>% arrange(padj) %>% head(50)
gsea_output1 <- res$combined |> as.data.frame()
gsea_output1$leadingEdge <- sapply(gsea_output1$leadingEdge, function(x) paste(x, collapse = ", "))
write.csv(gsea_output1, "edgeR_new_exploration/deg/Ribo_DMSO/GSEA/full_genes_list_GSEA.csv")
save(res, file = 'edgeR_new_exploration/deg/Ribo_DMSO/GSEA/GSEA.RData')

#' 
## ----Testing GSEA viz----------------------------------------------------------------------------------------------
# Function to barplot up and down GSEA results separatedly
plot_bulk_gsea_barplot_by_cluster <- function(gsea_df,cluster_col = "group",output_dir = "deg/GSEA/",label = "sCRNAseq") {
  library(ggpubr)
  library(ggplot2)
  # Create output dir
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  clusts <- sort(unique(gsea_df[[cluster_col]]))
  for (cl in clusts) {
    df_cl <- gsea_df[gsea_df[[cluster_col]] == cl, , drop = FALSE]
    # Preprocess for plotting (same as your function)
    tpm <- df_cl %>%mutate(adjPvalue = ifelse(padj <= 0.05, "significant", "non-significant"),
             pathway = gsub("HALLMARK_", "", pathway))  # Clean pathway names if Hallmark
    cols <- c("non-significant" = "black", "significant" = "#CD4071FF")
    p <- ggbarplot(tpm %>% filter(ES > 0),
                   x = "pathway", y = "NES",fill = "adjPvalue",color = "white",
                   palette = cols,sort.val = "asc",sort.by.groups = FALSE,
                   x.text.angle = 90,ylab = "NES",legend.title = "FDR < 0.05",rotate = TRUE,ggtheme = theme_minimal()) +
      theme(
        axis.text.x = element_text(size = 18, color = "black"),
        axis.title.x = element_text(size = 16, color = "black", face = "bold"),
        axis.text.y = element_text(size = 14, color = "black"),
        axis.title.y = element_text(size = 16, color = "black", face = "bold"),
        legend.text = element_text(size = 12),legend.title = element_text(size = 14))
    ggsave(
      filename = file.path(output_dir, paste0(label, "_cluster-", cl, "_Upregulated_GSEA.jpeg")),
      plot = p,width = 9,height = 9,dpi = 600,units = "in")}}

# Example of Vizualization
# res is the GSEA results from specific comparisons
# Re use this function for each result...
plot_bulk_gsea_barplot_by_cluster(
  gsea_df   = res$combined,   # or your per-cluster GSEA tibble
  cluster_col = "cluster",  output_dir = "edgeR_new_exploration/deg/Combo_DMSO/GSEA/barplots",label = "Combo_DMSO")

#' 
## ----GSEA visualization for publication----------------------------------------------------------------------------------
##--- Combo vs DMSO
load("edgeR_new_exploration/deg/Combo_DMSO/GSEA/GSEA.RData")
tpm_down <- res$combined %>%filter(padj <= 0.05, NES < 0) %>%
  mutate(pathway = gsub("HALLMARK_", "", pathway),pathway = gsub("_", " ", pathway),neg_log10_fdr = -log10(padj)) %>%
  arrange(NES)
tpm_down <- tpm_down %>%mutate(pathway = fct_reorder(pathway, NES))
# Viz
p <- ggplot(tpm_down, aes(x = NES, y = pathway)) +
  geom_col(aes(fill = neg_log10_fdr),width = 0.75,color = "white",linewidth = 0.3) +
  scale_fill_gradient(low = "#A6CEE3",high = "#08306B",name = expression(-log[10]("FDR"))) +
  geom_vline(xintercept = 0,linewidth = 0.6,color = "black") +
  labs(x = "NES",y = NULL,title = "BMS-986158/Ribociclib vs Vehicle") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 15, face = "bold", color = "black", hjust = 0.5),
    axis.text.y = element_text(size = 13, color = "black"),
    axis.text.x = element_text(size = 13, color = "black"),
    axis.title.x = element_text(size = 15, color = "black"),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 12),panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),plot.margin = margin(10, 20, 10, 10))
pdf("edgeR_new_exploration/deg/Combo_DMSO/GSEA/modified_down_barplot.pdf", width = 8, height = 6)
p
dev.off()

##--- BMS vs DMSO
load("edgeR_new_exploration/deg/BMS_DMSO/GSEA/GSEA.RData")
tpm_all <- res$combined %>% filter(padj <= 0.05) %>%mutate(pathway = gsub("HALLMARK_", "", pathway),pathway = gsub("_", " ", pathway),neg_log10_fdr = -log10(padj),direction = ifelse(NES > 0, "Up", "Down")) %>% arrange(NES) %>%mutate(pathway = forcats::fct_reorder(pathway, NES))
# Viz
p <- ggplot(tpm_all, aes(x = NES, y = pathway)) +
  geom_col(aes(fill = neg_log10_fdr),width = 0.75,color = "white",linewidth = 0.3) +
  scale_fill_gradient(low = "#EAF7E3",high = "#6FBF73",name = expression(-log[10]("FDR"))) +
  geom_vline(xintercept = 0, linewidth = 0.6, color = "black") +
  labs(x = "NES",y = NULL,title = "BMS-986158 vs Vehicle"
  ) + theme_minimal(base_size = 14)
pdf("edgeR_new_exploration/deg/BMS_DMSO/GSEA/modified_all_barplot.pdf", width = 8, height = 6)
p
dev.off()

#' 
#' # Volcano plots
## ----Volcano plots-------------------------------------------------------------------------------------------------
Combo <- read_xlsx("edgeR_new_exploration/deg/Combo_DMSO/CTX_FullTab_shrink.xlsx") # Combo vs DMSO
BMS <- read_xlsx("edgeR_new_exploration/deg/BMS_DMSO/CTX_FullTab_shrink.xlsx") # BMS vs DMSO
Ribo <- read_xlsx("edgeR_new_exploration/deg/Ribo_DMSO/CTX_FullTab_shrink.xlsx") # Ribo vs DMSO

# Temporary data frame for Combo vs DMSO
#--------------------------------
df <- Combo %>%
  mutate(LOG = -log10(padj), ABS = abs(log2FoldChange)) %>% mutate(Threshold = if_else(padj < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>% mutate(Direction = case_when(log2FoldChange > 0.3 & padj < 0.05 ~ "UpReg", log2FoldChange < -0.3 & padj < 0.05 ~ "DownReg"))

top_labelled <- df %>%group_by(Direction) %>%na.omit() %>%arrange(padj) %>%top_n(n = 10, wt = LOG) # top genes

## Volcano plot 
pdf("edgeR_new_exploration/deg/Combo_DMSO/Volcano_Plot_CTX.pdf",width=7,height=5,useDingbats=FALSE)
ggplot(df, aes(x = log2FoldChange, y = LOG)) +theme_classic(base_size = 13) +
  theme(axis.title.x = element_text(size = 20),axis.title.y = element_text(size = 20),
        axis.text.x  = element_text(size = 16),axis.text.y  = element_text(size = 16),
        axis.ticks   = element_line(linewidth = 0.8),axis.ticks.length = unit(0.25, "cm")) +
  geom_point(aes(color = Direction), size = 1, alpha = 0.8, show.legend = FALSE) +
  scale_color_manual(values = c(DownReg = "#D8CE7A", NS = "grey90", UpReg = "blue4")) +
  ggnewscale::new_scale_color() + geom_text_repel(data = top_labelled,aes(label = Gene, color = Direction),
    max.overlaps = Inf,size = 4, show.legend = FALSE) +
  geom_point(data = top_labelled,aes(fill = Direction),shape = 21,color = "black",stroke = 1,size = 2,alpha = 1,show.legend = FALSE) + 
  scale_fill_manual(values = c(DownReg = "#D8CE7A",UpReg = "blue4")) +
  scale_color_manual(values = c(DownReg = "black", UpReg = "black")) +
  geom_vline(xintercept = 0, linetype = 2, color = "red3") +
  geom_hline(yintercept = -log10(0.05), linetype = 2, color = "red3") +
  xlab("avg_log2FC") + ylab("-log10(adj p-value)") + scale_x_continuous(limits = c(-8, 8))
dev.off()

## Temporary dataframe for BMS vs DMSO
#--------------------------------
df <- BMS %>%mutate(LOG = -log10(padj), ABS = abs(log2FoldChange)) %>%
  mutate(Threshold = if_else(padj < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%mutate(Direction = case_when(log2FoldChange > 0.3 & padj < 0.05 ~ "UpReg", log2FoldChange < -0.3 & padj < 0.05 ~ "DownReg")) 

top_labelled <- df %>%group_by(Direction) %>%na.omit() %>%arrange(padj) %>%top_n(n = 10, wt = LOG)

## Volcano plot 
pdf("edgeR_new_exploration/deg/BMS_DMSO/Volcano_Plot_CTX.pdf",width=7,height=5,useDingbats=FALSE)
ggplot(df, aes(x = log2FoldChange, y = LOG)) +theme_classic(base_size = 13) +
  theme(axis.title.x = element_text(size = 20),axis.title.y = element_text(size = 20),
        axis.text.x  = element_text(size = 16),axis.text.y  = element_text(size = 16),
        axis.ticks   = element_line(linewidth = 0.8),axis.ticks.length = unit(0.25, "cm")) +
  geom_point(aes(color = Direction), size = 1, alpha = 0.8, show.legend = FALSE) +
  scale_color_manual(values = c(DownReg = "#D8CE7A", NS = "grey90", UpReg = "#A5D996")) +
  ggnewscale::new_scale_color() +geom_text_repel(data = top_labelled,aes(label = Gene, color = Direction),
    max.overlaps = Inf,size = 4, show.legend = FALSE) +
  geom_point(data = top_labelled,aes(fill = Direction),shape = 21,color = "black",stroke = 1,size = 2,alpha = 1,show.legend = FALSE) + 
  scale_fill_manual(values = c(DownReg = "#D8CE7A",UpReg = "#A5D996")) +
  scale_color_manual(values = c(DownReg = "black", UpReg = "black")) +
  geom_vline(xintercept = 0, linetype = 2, color = "red3") +
  geom_hline(yintercept = -log10(0.05), linetype = 2, color = "red3") +
  xlab("avg_log2FC") + ylab("-log10(adj p-value)") + scale_x_continuous(limits = c(-8, 8))
dev.off()

## Temporary dataframe for Ribo vs DMSO
#--------------------------------
df <- Ribo %>%mutate(LOG = -log10(padj), ABS = abs(log2FoldChange)) %>%
  mutate(Threshold = if_else(padj < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%mutate(Direction = case_when(log2FoldChange > 0.3 & padj < 0.05 ~ "UpReg", log2FoldChange < -0.3 & padj < 0.05 ~ "DownReg")) ### why 0.3

top_labelled <- df %>%group_by(Direction) %>%na.omit() %>%arrange(padj) %>%top_n(n = 10, wt = LOG)

## Volcano plot 
pdf("edgeR_new_exploration/deg/Ribo_DMSO/Volcano_Plot_CTX.pdf",width=7,height=5,useDingbats=FALSE)
ggplot(df, aes(x = log2FoldChange, y = LOG)) +theme_classic(base_size = 13) +
  theme( axis.title.x = element_text(size = 20),axis.title.y = element_text(size = 20),
    axis.text.x  = element_text(size = 16),axis.text.y  = element_text(size = 16),axis.ticks   = element_line(linewidth = 0.8),axis.ticks.length = unit(0.25, "cm")) +
  geom_point(aes(color = Direction), size = 1, alpha = 0.8, show.legend = FALSE) +
  scale_color_manual(values = c(DownReg = "#D8CE7A", NS = "grey90", UpReg = "#8259D9")) +
  ggnewscale::new_scale_color() +
  geom_text_repel(data = top_labelled,aes(label = Gene, color = Direction),max.overlaps = Inf,size = 4, show.legend = FALSE) +
  geom_point(data = top_labelled,aes(fill = Direction),shape = 21,color = "black",stroke = 1,size = 2,alpha = 1,show.legend = FALSE) + 
  scale_fill_manual(values = c(DownReg = "#D8CE7A",UpReg = "#8259D9")) +
  scale_color_manual(values = c(DownReg = "black", UpReg = "black")) +
  geom_vline(xintercept = 0, linetype = 2, color = "red3") +
  geom_hline(yintercept = -log10(0.05), linetype = 2, color = "red3") +
  xlab("avg_log2FC") + ylab("-log10(adj p-value)") + scale_x_continuous(limits = c(-8, 8))
dev.off()

#' 
#' # Heatmaps
## ----heatmpas on Combo vs DMSO-------------------------------------------------------------------------------------
# pheatmap of Combo vs DMSO using vst normalized scores 
#Combo_DMSO <- read_xlsx("edgeR_new_exploration/deg/Combo_DMSO/CTX_DGE_shrink.xlsx")
load("edgeR_new_exploration/deg/dds.RData")
vsd <- vst(dds, blind = FALSE)
mat <- assay(vsd)
# Filtering Combo and DMSO samples
mat <- mat[rownames(mat) %in% Combo_DMSO$Gene, ]
mat <- mat[, c(1:3, 10:12)]

# loading metadata
anno<- read.csv("edgeR_new_exploration/metadata.csv")
anno$Subject_ID<- NULL
anno$Sex<- NULL
anno<- anno|> filter(Condition %in% c("DMSO", "Combo"))
anno$Condition <- factor(anno$Condition, levels = c("DMSO","Combo"))

# Viz of testing heatmap
my_colour = list(Condition = c(DMSO = "#D8CE7A", Combo = "blue4"))
pdf("edgeR_new_exploration/deg/Combo_BMS/DEGs_Heatmap_CTX.pdf",width=6,height=7,useDingbats=FALSE)
pheatmap::pheatmap(mat,scale="row",
         show_rownames = F,annotation= anno,annotation_colors = my_colour,
         colorRampPalette(rev(viridis(100)))(100),cellwidth = 12,
         #cutree_cols = 2,
         angle_col = 45,cluster_cols = F)
dev.off()

# Heatmap of Combo vs DMSO showing AMPK significant genes
# Scaling mat for viz
mat_scaled <- t(scale(t(mat)))
mat_scaled[is.na(mat_scaled)] <- 0

ha <- HeatmapAnnotation(
  Comparison = colnames(mat_scaled),
  col = list(Comparison = c(
      DMSO_1 = "#D8CE7A",DMSO_2 = "#D8CE7A",DMSO_3 = "#D8CE7A",
      Combo_1 = "blue4",Combo_3 = "blue4",Combo_4 = "blue4")),show_legend = T)

genes_to_show <- c("Aak1", "Araf", " Bcr", " Cdc42bpa", " Mtor", "Dclk1","Ern1", "Gsk3b", "Lmtk2", "Mtor", "Pikfyve", 
                   "Prkaa1","Prkag1", "Prkx", "Pskh1", "Ripk2", "Rps6ka2", "Sik2","Srpk2", "Stk36", "Taok2", 
                   "Trio", "Ttbk1", 'Ttbk2',"Wnk3", 'Wnk4') # AMPK DEGs

lab <-which(rownames(mat_scaled) %in% genes_to_show)
hr = rowAnnotation(foo = anno_mark(at = which(rownames(mat_scaled) %in% genes_to_show),
                                    labels = rownames(mat_scaled)[rownames(mat_scaled)%in%genes_to_show],
                                   labels_gp = gpar(fontsize = 10)))

# First viz
pdf("edgeR_new_exploration/deg/Combo_DMSO/DGE_heatmap_foldCs.pdf", width=6, height=8)
Heatmap(mat_scaled, cluster_rows=T,show_row_dend=T,cluster_columns=T,col=circlize::colorRamp2(c(-2,0,2),c("cyan", "black", "yellow")), show_column_names=F, show_row_names=F, column_title=NULL,top_annotation=ha, 
        right_annotation = hr,row_km = 4,row_names_max_width = unit(10,"in"),row_title_gp = gpar(fontsize = 15),
        column_title_gp = gpar(fontsize = 20),column_names_gp = gpar(fontsize =10),row_names_gp = gpar(fontsize =3)) 
dev.off()

# Second viz
pdf("edgeR_new_exploration/deg/Combo_DMSO/DGE_heatmap_foldCs_1.pdf", width=6, height=8)
Heatmap(mat_scaled, cluster_rows=T,show_row_dend=T,cluster_columns=T,col=circlize::colorRamp2(c(-2,0,2),c("#053061", "#FDF5E6", "#67001F")),show_column_names=F, show_row_names=F, column_title=NULL,top_annotation=ha, 
        right_annotation = hr,row_km = 4,row_names_max_width = unit(10,"in"),row_title_gp = gpar(fontsize = 15),
        column_title_gp = gpar(fontsize = 20),column_names_gp = gpar(fontsize =10),row_names_gp = gpar(fontsize =3)) 
dev.off()

# Third viz 
pdf("edgeR_new_exploration/deg/Combo_DMSO/DGE_heatmap_foldCs_2.pdf", width=6, height=8)
Heatmap(mat_scaled, cluster_rows=T,show_row_dend=T,cluster_columns=T,col=circlize::colorRamp2(c(-2, -1, 0, 1, 2),viridis(5)),show_column_names=F, show_row_names=F, 
        column_title=NULL,top_annotation=ha, right_annotation = hr,
        row_km = 4,  row_names_max_width = unit(10,"in"),row_title_gp = gpar(fontsize = 15),
        column_title_gp = gpar(fontsize = 20),column_names_gp = gpar(fontsize =10),row_names_gp = gpar(fontsize =3)) 
dev.off()

# Fourth viz
pdf("edgeR_new_exploration/deg/Combo_DMSO/DGE_heatmap_foldCs_3.pdf", width=6, height=8)
Heatmap(mat_scaled, cluster_rows=T,show_row_dend=T,cluster_columns=T,col=circlize::colorRamp2(c(-2, -1, 0, 1, 2),c("#2166AC", "#67A9CF", "#F7F7F7", "#EF8A62", "#B2182B")),show_column_names=F, 
        show_row_names=F, column_title=NULL,top_annotation=ha, right_annotation = hr,
        row_km = 4,row_names_max_width = unit(10,"in"),row_title_gp = gpar(fontsize = 15),
        column_title_gp = gpar(fontsize = 20),column_names_gp = gpar(fontsize =10),row_names_gp = gpar(fontsize =3)) 
dev.off()

#' 
## ----GSEA GO ontology barplots-------------------------------------------------------------------------------------
Combo<- read_xlsx("edgeR_new_exploration/deg/Combo_DMSO/CTX_FullTab_shrink.xlsx") # Combo vs DMSO
Ribo<- read_xlsx("edgeR_new_exploration/deg/Ribo_DMSO/CTX_FullTab_shrink.xlsx") # Ribo vs DMSO
BMS<- read_xlsx("edgeR_new_exploration/deg/BMS_DMSO/CTX_FullTab_shrink.xlsx") # BMS vs DMSO

#-----------------------------
# GSEA-GO on Combo vs DMSO
logFC <- Combo$log2FoldChange
names(logFC) <- Combo$Gene
# Map SYMBOL -> ENTREZID (mouse)
conv <- bitr(names(logFC),fromType = "SYMBOL",toType   = "ENTREZID",OrgDb = org.Mm.eg.db)
df <- data.frame(SYMBOL = names(logFC), logFC = as.numeric(logFC))
df <- dplyr::inner_join(df, conv, by = "SYMBOL")
df <- df %>% dplyr::group_by(ENTREZID) %>%dplyr::summarise(logFC = logFC[which.max(abs(logFC))], .groups = "drop")
geneList <- df$logFC
names(geneList) <- df$ENTREZID
geneList <- sort(geneList, decreasing = TRUE)

#### Fixing gene list
#geneList <- geneList[is.finite(geneList)]
#geneList <- geneList[!is.na(names(geneList))]
#geneList <- geneList[!duplicated(names(geneList))]

set.seed(708)
## Running GSEA GO(BP): Change to MF when needed
g <- gseGO(geneList = geneList,OrgDb= org.Mm.eg.db, keyType= "ENTREZID", ont= "BP", minGSSize=10, maxGSSize=500, pvalueCutoff = 1,verbose= FALSE, nPerm = 10000, eps= 0)
if (is.null(g) || nrow(g@result) == 0) return(NULL)
out <- g@result %>%dplyr::select(ID, Description, setSize, NES, pvalue, qvalue, core_enrichment) %>% dplyr::filter(qvalue < 0.05)
out
## viz 
top5_up <- out %>%filter(NES > 0) %>%slice_max(order_by = NES, n = 20, with_ties = FALSE) %>% ungroup()
top5_down <- out %>%filter(NES < 0) %>%slice_min(order_by = NES, n = 20, with_ties = FALSE) %>%ungroup()
top5_both <- bind_rows(top5_up, top5_down) %>%mutate(Direction = ifelse(NES > 0, "Up", "Down"))
top5_both <- bind_rows(top5_up, top5_down) %>%distinct(ID, .keep_all = TRUE) %>%arrange(NES) %>%                      mutate(Description = forcats::fct_inorder(Description))

pad <- 0.03 * diff(range(top5_both$NES, na.rm = TRUE))
p<- ggplot(top5_both, aes(x = NES, y = Description, fill = NES > 0)) + geom_col(width = 0.75) + geom_vline(xintercept = 0, linewidth = 0.4) +
  # GO IDs inside bars
  geom_text(aes(label = ID), hjust = ifelse(top5_both$NES > 0, 1.05, -0.05), color = "white", size = 4) +
  geom_text(aes(x = ifelse(NES > 0, -pad, pad),label = Description),hjust = ifelse(top5_both$NES > 0, 1, 0),size = 3) +
  tidytext::scale_y_reordered() + scale_fill_manual(values = c(`TRUE` = "firebrick", `FALSE` = "firebrick3")) + ##4575B4
  theme_classic(base_size = 12) + theme(plot.title = element_text(hjust = 0.5, size = 16),legend.position = "none",
    axis.line.y  = element_blank(),axis.ticks.y = element_blank(),axis.text.y  = element_blank(),
    axis.text.x  = element_text(size = 14),axis.title.x = element_text(size = 16),
    plot.margin = margin(5.5, 80, 5.5, 80)
  ) + coord_cartesian(clip = "off") +labs(x = "NES", y = NULL,title = "GSEA (BP) | Combo vs DMSO")

pdf("edgeR_new_exploration/deg/Combo_DMSO/GO_GSEA/BP_GO_GSEA_plot.pdf", width = 9, height = 7)
p
dev.off()

#-----------------------------
# GSEA-GO on BMS vs DMSO
logFC <- BMS$log2FoldChange
names(logFC) <- BMS$Gene
# Map SYMBOL -> ENTREZID (mouse)
conv <- bitr(names(logFC),fromType = "SYMBOL",toType   = "ENTREZID",OrgDb = org.Mm.eg.db)
df <- data.frame(SYMBOL = names(logFC), logFC = as.numeric(logFC))
df <- dplyr::inner_join(df, conv, by = "SYMBOL")
df <- df %>% dplyr::group_by(ENTREZID) %>%dplyr::summarise(logFC = logFC[which.max(abs(logFC))], .groups = "drop")
geneList <- df$logFC
names(geneList) <- df$ENTREZID
geneList <- sort(geneList, decreasing = TRUE)

set.seed(708)
## Running GSEA GO(BP): Change to MF when needed
g <- gseGO(geneList = geneList,OrgDb= org.Mm.eg.db, keyType= "ENTREZID", ont= "BP", minGSSize=10, maxGSSize=500, pvalueCutoff = 1,verbose= FALSE, nPerm = 10000, eps= 0)
if (is.null(g) || nrow(g@result) == 0) return(NULL)
out <- g@result %>%dplyr::select(ID, Description, setSize, NES, pvalue, qvalue, core_enrichment) %>% dplyr::filter(qvalue < 0.05)
out
## viz 
top5_up <- out %>%filter(NES > 0) %>%slice_max(order_by = NES, n = 20, with_ties = FALSE) %>% ungroup()
top5_down <- out %>%filter(NES < 0) %>%slice_min(order_by = NES, n = 20, with_ties = FALSE) %>%ungroup()
top5_both <- bind_rows(top5_up, top5_down) %>%mutate(Direction = ifelse(NES > 0, "Up", "Down"))
top5_both <- bind_rows(top5_up, top5_down) %>%distinct(ID, .keep_all = TRUE) %>%arrange(NES) %>%                      mutate(Description = forcats::fct_inorder(Description))

pad <- 0.03 * diff(range(top5_both$NES, na.rm = TRUE))
p<- ggplot(top5_both, aes(x = NES, y = Description, fill = NES > 0)) + geom_col(width = 0.75) + geom_vline(xintercept = 0, linewidth = 0.4) +
  # GO IDs inside bars
  geom_text(aes(label = ID), hjust = ifelse(top5_both$NES > 0, 1.05, -0.05), color = "white", size = 4) +
  geom_text(aes(x = ifelse(NES > 0, -pad, pad),label = Description),hjust = ifelse(top5_both$NES > 0, 1, 0),size = 3) +
  tidytext::scale_y_reordered() + scale_fill_manual(values = c(`TRUE` = "firebrick", `FALSE` = "firebrick3")) + ##4575B4
  theme_classic(base_size = 12) + theme(plot.title = element_text(hjust = 0.5, size = 16),legend.position = "none",
    axis.line.y  = element_blank(),axis.ticks.y = element_blank(),axis.text.y  = element_blank(),
    axis.text.x  = element_text(size = 14),axis.title.x = element_text(size = 16),
    plot.margin = margin(5.5, 80, 5.5, 80)
  ) + coord_cartesian(clip = "off") +labs(x = "NES", y = NULL,title = "GSEA (BP) | BMS vs DMSO")

pdf("edgeR_new_exploration/deg/BMS_DMSO/GO_GSEA/BP_GO_GSEA_plot.pdf", width = 9, height = 7)
p
dev.off()

#-----------------------------
# GSEA-GO on Ribo vs DMSO
logFC <- Ribo$log2FoldChange
names(logFC) <- Ribo$Gene
# Map SYMBOL -> ENTREZID (mouse)
conv <- bitr(names(logFC),fromType = "SYMBOL",toType   = "ENTREZID",OrgDb = org.Mm.eg.db)
df <- data.frame(SYMBOL = names(logFC), logFC = as.numeric(logFC))
df <- dplyr::inner_join(df, conv, by = "SYMBOL")
df <- df %>% dplyr::group_by(ENTREZID) %>%dplyr::summarise(logFC = logFC[which.max(abs(logFC))], .groups = "drop")
geneList <- df$logFC
names(geneList) <- df$ENTREZID
geneList <- sort(geneList, decreasing = TRUE)

set.seed(708)
## Running GSEA GO(BP): Change to MF when needed
g <- gseGO(geneList = geneList,OrgDb= org.Mm.eg.db, keyType= "ENTREZID", ont= "BP", minGSSize=10, maxGSSize=500, pvalueCutoff = 1,verbose= FALSE, nPerm = 10000, eps= 0)
if (is.null(g) || nrow(g@result) == 0) return(NULL)
out <- g@result %>%dplyr::select(ID, Description, setSize, NES, pvalue, qvalue, core_enrichment) %>% dplyr::filter(qvalue < 0.05)
out
## viz 
top5_up <- out %>%filter(NES > 0) %>%slice_max(order_by = NES, n = 20, with_ties = FALSE) %>% ungroup()
top5_down <- out %>%filter(NES < 0) %>%slice_min(order_by = NES, n = 20, with_ties = FALSE) %>%ungroup()
top5_both <- bind_rows(top5_up, top5_down) %>%mutate(Direction = ifelse(NES > 0, "Up", "Down"))
top5_both <- bind_rows(top5_up, top5_down) %>%distinct(ID, .keep_all = TRUE) %>%arrange(NES) %>%                      mutate(Description = forcats::fct_inorder(Description))

pad <- 0.03 * diff(range(top5_both$NES, na.rm = TRUE))
p<- ggplot(top5_both, aes(x = NES, y = Description, fill = NES > 0)) + geom_col(width = 0.75) + geom_vline(xintercept = 0, linewidth = 0.4) +
  # GO IDs inside bars
  geom_text(aes(label = ID), hjust = ifelse(top5_both$NES > 0, 1.05, -0.05), color = "white", size = 4) +
  geom_text(aes(x = ifelse(NES > 0, -pad, pad),label = Description),hjust = ifelse(top5_both$NES > 0, 1, 0),size = 3) +
  tidytext::scale_y_reordered() + scale_fill_manual(values = c(`TRUE` = "firebrick", `FALSE` = "firebrick3")) + ##4575B4
  theme_classic(base_size = 12) + theme(plot.title = element_text(hjust = 0.5, size = 16),legend.position = "none",
    axis.line.y  = element_blank(),axis.ticks.y = element_blank(),axis.text.y  = element_blank(),
    axis.text.x  = element_text(size = 14),axis.title.x = element_text(size = 16),
    plot.margin = margin(5.5, 80, 5.5, 80)
  ) + coord_cartesian(clip = "off") +labs(x = "NES", y = NULL,title = "GSEA (BP) | Ribo vs DMSO")

pdf("edgeR_new_exploration/deg/Ribo_DMSO/GO_GSEA/BP_GO_GSEA_plot.pdf", width = 9, height = 7)
p
dev.off()

#' 
## ----r session-info------------------------------------------------------------------------------------------------
sessioninfo::session_info()


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
# [1] grid      stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#  [1] org.Mm.eg.db_3.22.0         AnnotationDbi_1.72.0        sessioninfo_1.2.3           ComplexHeatmap_2.26.0      
#  [5] readxl_1.4.5                clusterProfiler_4.18.4      fgsea_1.36.2                scToppR_0.99.4             
#  [9] writexl_1.5.4               tidyplots_0.4.0             janitor_2.2.1               broom_1.0.12               
# [13] emmeans_2.0.3               limma_3.66.0                viridis_0.6.5               viridisLite_0.4.3          
# [17] sva_3.56.0                  BiocParallel_1.44.0         genefilter_1.90.0           mgcv_1.9-3                 
# [21] nlme_3.1-168                pheatmap_1.0.13             DESeq2_1.48.2               SummarizedExperiment_1.38.1
# [25] Biobase_2.68.0              MatrixGenerics_1.20.0       matrixStats_1.5.0           GenomicRanges_1.60.0       
# [29] GenomeInfoDb_1.44.3         IRanges_2.44.0              S4Vectors_0.48.0            BiocGenerics_0.54.1        
# [33] generics_0.1.4              future.apply_1.20.2         future_1.70.0               preprocessCore_1.71.2      
# [37] lubridate_1.9.5             forcats_1.0.1               stringr_1.6.0               dplyr_1.2.1                
# [41] purrr_1.2.2                 readr_2.2.0                 tidyr_1.3.2                 tibble_3.3.1               
# [45] tidyverse_2.0.0             RColorBrewer_1.1-3          data.table_1.18.2.1         ggpubr_0.6.3               
# [49] ggrepel_0.9.8               ggplot2_4.0.2              
# 
# loaded via a namespace (and not attached):
#   [1] splines_4.5.1           ggplotify_0.1.3         cellranger_1.1.0        R.oo_1.27.1            
#   [5] polyclip_1.10-7         XML_3.99-0.23           lifecycle_1.0.5         httr2_1.2.2            
#   [9] rstatix_0.7.3           doParallel_1.0.17       edgeR_4.8.2             MASS_7.3-65            
#  [13] globals_0.19.1          lattice_0.22-7          backports_1.5.1         magrittr_2.0.4         
#  [17] openxlsx_4.2.8.1        otel_0.2.0              ggtangle_0.1.1          zip_2.3.3              
#  [21] cowplot_1.2.0           DBI_1.3.0               abind_1.4-8             R.utils_2.13.0         
#  [25] yulab.utils_0.2.4       tweenr_2.0.3            rappdirs_0.3.4          circlize_0.4.18        
#  [29] gdtools_0.5.0           GenomeInfoDbData_1.2.14 enrichplot_1.30.4       listenv_0.10.1         
#  [33] tidytree_0.4.7          annotate_1.88.0         parallelly_1.46.1       codetools_0.2-20       
#  [37] DelayedArray_0.36.0     ggforce_0.5.0           DOSE_4.4.0              shape_1.4.6.1          
#  [41] tidyselect_1.2.1        aplot_0.2.9             UCSC.utils_1.4.0        farver_2.1.2           
#  [45] Seqinfo_1.0.0           jsonlite_2.0.0          GetoptLong_1.1.1        Formula_1.2-5          
#  [49] iterators_1.0.14        survival_3.8-3          systemfonts_1.3.2       foreach_1.5.2          
#  [53] tools_4.5.1             ggnewscale_0.5.2        treeio_1.34.0           Rcpp_1.1.1             
#  [57] glue_1.8.0              gridExtra_2.3           SparseArray_1.10.8      xfun_0.55              
#  [61] qvalue_2.42.0           withr_3.0.2             fastmap_1.2.0           digest_0.6.39          
#  [65] timechange_0.4.0        R6_2.6.1                gridGraphics_0.5-1      estimability_1.5.1     
#  [69] colorspace_2.1-2        GO.db_3.22.0            dichromat_2.0-0.1       RSQLite_2.4.6          
#  [73] R.methodsS3_1.8.2       fontLiberation_0.1.0    htmlwidgets_1.6.4       httr_1.4.8             
#  [77] S4Arrays_1.10.1         scatterpie_0.2.6        pkgconfig_2.0.3         gtable_0.3.6           
#  [81] blob_1.3.0              S7_0.2.1                XVector_0.50.0          htmltools_0.5.9        
#  [85] fontBitstreamVera_0.1.1 carData_3.0-6           clue_0.3-68             scales_1.4.0           
#  [89] png_0.1-9               snakecase_0.11.1        ggfun_0.2.0             knitr_1.51             
#  [93] rstudioapi_0.18.0       rjson_0.2.23            tzdb_0.5.0              reshape2_1.4.5         
#  [97] coda_0.19-4.1           GlobalOptions_0.1.4     cachem_1.1.0            parallel_4.5.1         
# [101] pillar_1.11.1           vctrs_0.7.2             tidydr_0.0.6            car_3.1-5              
# [105] cluster_2.1.8.2         xtable_1.8-8            evaluate_1.0.5          mvtnorm_1.3-6          
# [109] cli_3.6.6               locfit_1.5-9.12         compiler_4.5.1          rlang_1.1.7            
# [113] crayon_1.5.3            ggsignif_0.6.4          plyr_1.8.9              fs_2.0.1               
# [117] ggiraph_0.9.6           stringi_1.8.7           Biostrings_2.78.0       lazyeval_0.2.3         
# [121] fontquiver_0.2.1        GOSemSim_2.36.0         Matrix_1.7-3            hms_1.1.4              
# [125] patchwork_1.3.2         bit64_4.6.0-1           KEGGREST_1.50.0         statmod_1.5.1          
# [129] igraph_2.2.3            memoise_2.0.1           ggtree_4.0.4            fastmatch_1.1-8        
# [133] bit_4.6.0               gson_0.1.0              ape_5.8-1              
