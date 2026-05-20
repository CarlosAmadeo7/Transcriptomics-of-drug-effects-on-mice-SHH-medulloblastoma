#' ---
#' title: "Network enrichment analysis"
#' output: html_document
#' date: "2025-11-25"
#' ---
#' 
## ----setup, include=FALSE------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

#' 
## ----Libraries-----------------------------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(readxl);library(dplyr);library(clusterProfiler);library(org.Mm.eg.db);
  library(ReactomePA);library(enrichplot);library(ggplot2);library(DOSE);library(sessioninfo)
})

#' 
## ----Loading files-------------------------------------------------------------------------------------------------
Combo <- read_xlsx("edgeR_new_exploration/deg/Combo_DMSO/CTX_DGE_shrink.xlsx") %>% 
  mutate(Condition = "BMS-986158/Ribociclib")
BMS <- read_xlsx("edgeR_new_exploration/deg/BMS_DMSO/CTX_DGE_shrink.xlsx") %>% 
  mutate(Condition = "BMS-986158")
Ribo <- read_xlsx("edgeR_new_exploration/deg/Ribo_DMSO/CTX_DGE_shrink.xlsx") %>%
  mutate(Condition = "Ribociclib")

df <- bind_rows(Combo, BMS, Ribo) %>%
  mutate(LOG = -log10(padj), ABS = abs(log2FoldChange)) %>%
  mutate(Threshold = if_else(padj < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
  mutate(Direction = case_when(log2FoldChange > 0.3 & padj < 0.05 ~ "UpReg", log2FoldChange < -0.3 & padj < 0.05 ~ "DownReg"))

#' 
## ----Enrichment plots of Upregulated genes-------------------------------------------------------------------------
## Enrichment plot of Upregulated genes 
df<- df |> filter(Direction == "UpReg")
markers_list <- split(df$Gene, df$Condition)

for(i in 1:length(markers_list)){
  markers_list[[i]] <- bitr(markers_list[[i]], fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")}

genelist <- list()
for(i in 1:length(markers_list)){
  genelist[[i]] <- markers_list[[i]]$ENTREZID}

names(genelist) <- names(markers_list)

# GO_BP
markers_GOBP <- compareCluster(geneClusters = genelist, fun = "enrichGO", ont="BP", OrgDb="org.Mm.eg.db")
# GO_MF
markers_GOMF <- compareCluster(geneClusters = genelist, fun = "enrichGO", ont="MF", OrgDb="org.Mm.eg.db")
# KEGG
markers_kegg <- compareCluster(geneClusters = genelist, fun="enrichKEGG", organism="mmu", keyType="kegg")
# WikiPathway
markers_WP <- compareCluster(geneClusters = genelist, fun="enrichWP", organism="Mus musculus")
# Reactomeでclusterprofiler
markers_Reactome <- compareCluster(genelist, fun = "enrichPathway", organism = "mouse")

dotplot(markers_GOBP, showCategory=7) 
dotplot(markers_GOMF, showCategory=7) 
dotplot(markers_kegg, showCategory=7) 
dotplot(markers_Reactome, showCategory=7) 
dotplot(markers_WP, showCategory=7) 

save(markers_GOBP, file = "edgeR_new_exploration/deg/Enrichmnet_plots/Up_markers_GOBP.RData")
save(markers_GOMF, file = "edgeR_new_exploration/deg/Enrichmnet_plots/Up_markers_GOMF.RData")
save(markers_kegg, file ="edgeR_new_exploration/deg/Enrichmnet_plots/Up_markers_kegg.RData")
save(markers_WP, file ="edgeR_new_exploration/deg/Enrichmnet_plots/Up_markers_WP.RData")
save(markers_Reactome, file ="edgeR_new_exploration/deg/Enrichmnet_plots/Up_markers_Reactome.RData")

## Viz
ct_cols <- c("BMS-986158/Ribociclib" = "royalblue4","BMS-986158" = "#A5D996","Ribociclib" = "#8259D9")
# GO_BP
p<- enrichplot::cnetplot(markers_GOBP,showCategory = 6,size_category = 3, node_label= 'category', size_edge = 0.01, size_item = 1)

pdf("edgeR_new_exploration/deg/Enrichmnet_plots/Up_markers_GOBP.pdf",height = 9, width = 10)
p + scale_fill_manual(values = ct_cols) + theme(text = element_text(size = 10),legend.text = element_text(size = 10),legend.title = element_blank())
dev.off()

# GO_MF
p<- enrichplot::cnetplot(markers_GOMF,showCategory = 7,size_category = 3,node_label= 'category', size_edge = 0.01, size_item = 1)
pdf("edgeR_new_exploration/deg/Enrichmnet_plots/Up_markers_GOMF.pdf",height = 9, width = 10)
p + scale_fill_manual(values = ct_cols) + theme(text = element_text(size = 10),legend.text = element_text(size = 10),legend.title = element_blank())
dev.off()

# KEEG
p<- enrichplot::cnetplot(markers_kegg,showCategory = 6,size_category = 3,node_label= 'category', size_edge = 0.01, size_item = 1)
pdf("edgeR_new_exploration/deg/Enrichmnet_plots/Up_markers_KEEG.pdf",height = 9, width = 10)
p + scale_fill_manual(values = ct_cols) + theme(text = element_text(size = 10),legend.text = element_text(size = 10),legend.title = element_blank())
dev.off()

# REACTOME
p<- enrichplot::cnetplot(markers_Reactome, showCategory = 7, size_category = 3,node_label= 'category', size_edge = 0.01, size_item = 1)
pdf("edgeR_new_exploration/deg/Enrichmnet_plots/Up_markers_REACTOME.pdf",height = 9, width = 10)
p + scale_fill_manual(values = ct_cols) + theme(text = element_text(size = 10), legend.text = element_text(size = 10),legend.title = element_blank())
dev.off()

# WikiPathways
p<- enrichplot::cnetplot(markers_WP,showCategory = 7,size_category = 3,node_label= 'category', size_edge = 0.01, size_item = 1)
pdf("edgeR_new_exploration/deg/Enrichmnet_plots/Up_markers_WP.pdf",height = 9, width = 10)
p + scale_fill_manual(values = ct_cols) + theme(text = element_text(size = 10),legend.text = element_text(size = 10),legend.title = element_blank())
dev.off()

#' 
## ----Enrichment plots of DownRegulated genes-----------------------------------------------------------------------
df <- bind_rows(Combo, BMS, Ribo) %>% mutate(LOG = -log10(padj), ABS = abs(log2FoldChange)) %>%
  mutate(Threshold = if_else(padj < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
  mutate(Direction = case_when(log2FoldChange > 0.3 & padj < 0.05 ~ "UpReg", log2FoldChange < -0.3 & padj < 0.05 ~ "DownReg"))

df<- df |> filter(Direction == "DownReg")
markers_list <- split(df$Gene, df$Condition)
## Enrichment analysis
for(i in 1:length(markers_list)){
  markers_list[[i]] <- bitr(markers_list[[i]], fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")}
genelist <- list()
for(i in 1:length(markers_list)){
  genelist[[i]] <- markers_list[[i]]$ENTREZID}
names(genelist) <- names(markers_list)

# GO_BP
markers_GOBP <- compareCluster(geneClusters = genelist, fun = "enrichGO", ont="BP", OrgDb="org.Mm.eg.db")
# GO_MF
markers_GOMF <- compareCluster(geneClusters = genelist, fun = "enrichGO", ont="MF", OrgDb="org.Mm.eg.db")
# KEGG
markers_kegg <- compareCluster(geneClusters = genelist, fun="enrichKEGG", organism="mmu", keyType="kegg")
# WikiPathway
markers_WP <- compareCluster(geneClusters = genelist, fun="enrichWP", organism="Mus musculus")
# Reactomeでclusterprofiler
markers_Reactome <- compareCluster(genelist, fun = "enrichPathway", organism = "mouse")

dotplot(markers_GOBP, showCategory=7) 
dotplot(markers_GOMF, showCategory=7) 
dotplot(markers_kegg, showCategory=7) 
dotplot(markers_Reactome, showCategory=7) 
dotplot(markers_WP, showCategory=7) 

save(markers_GOBP, file = "edgeR_new_exploration/deg/Enrichmnet_plots/Down_markers_GOBP.RData")
save(markers_GOMF, file = "edgeR_new_exploration/deg/Enrichmnet_plots/Down_markers_GOMF.RData")
save(markers_kegg, file ="edgeR_new_exploration/deg/Enrichmnet_plots/Down_markers_kegg.RData")
save(markers_WP, file ="edgeR_new_exploration/deg/Enrichmnet_plots/Down_markers_WP.RData")
save(markers_Reactome, file ="edgeR_new_exploration/deg/Enrichmnet_plots/Down_markers_Reactome.RData")

## Viz
ct_cols <- c("BMS-986158/Ribociclib" = "royalblue4","BMS-986158" = "#A5D996","Ribociclib" = "#8259D9")

# GO_BP
p<- enrichplot::cnetplot(markers_GOBP,showCategory = 5,size_category = 3,node_label= 'category', size_edge = 0.01, size_item = 1)
pdf("edgeR_new_exploration/deg/Enrichmnet_plots/Down_markers_GOBP.pdf",height = 9, width = 10)
p + scale_fill_manual(values = ct_cols) + theme(text = element_text(size = 10), legend.text = element_text(size = 10),legend.title = element_blank())
dev.off()

# GO_MF
p<- enrichplot::cnetplot(markers_GOMF,showCategory = 7,size_category = 3,node_label= 'category', size_edge = 0.01, size_item = 1)
pdf("edgeR_new_exploration/deg/Enrichmnet_plots/Down_markers_GOMF.pdf",height = 9, width = 10)
p + scale_fill_manual(values = ct_cols) + theme(text = element_text(size = 10), legend.text = element_text(size = 10),legend.title = element_blank())
dev.off()

# KEEG
p<- enrichplot::cnetplot(markers_kegg,showCategory = 6,size_category = 3,node_label= 'category', size_edge = 0.01, size_item = 1)
pdf("edgeR_new_exploration/deg/Enrichmnet_plots/Down_markers_KEEG.pdf",height = 9, width = 10)
p + scale_fill_manual(values = ct_cols) + theme(text = element_text(size = 10),legend.text = element_text(size = 10),legend.title = element_blank())
dev.off()

# REACTOME
p<- enrichplot::cnetplot(markers_Reactome,showCategory = 7,size_category = 3,node_label= 'category', size_edge = 0.01, size_item = 1)
pdf("edgeR_new_exploration/deg/Enrichmnet_plots/Down_markers_REACTOME.pdf",height = 9, width = 10)
p + scale_fill_manual(values = ct_cols) + theme(text = element_text(size = 10), legend.text = element_text(size = 10),legend.title = element_blank())
dev.off()

# WikiPathways
p<- enrichplot::cnetplot(markers_WP,showCategory = 7,size_category = 3,node_label= 'category', size_edge = 0.01, size_item = 1)
pdf("edgeR_new_exploration/deg/Enrichmnet_plots/Down_markers_WP.pdf",height = 9, width = 10)
p + scale_fill_manual(values = ct_cols) + theme(text = element_text(size = 10), legend.text = element_text(size = 10),legend.title = element_blank())
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
# [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#  [1] sessioninfo_1.2.3      DOSE_4.4.0             ggplot2_4.0.2          enrichplot_1.30.4     
#  [5] ReactomePA_1.52.0      org.Mm.eg.db_3.22.0    AnnotationDbi_1.72.0   IRanges_2.44.0        
#  [9] S4Vectors_0.48.0       Biobase_2.68.0         BiocGenerics_0.54.1    generics_0.1.4        
# [13] clusterProfiler_4.18.4 dplyr_1.2.1            readxl_1.4.5          
# 
# loaded via a namespace (and not attached):
#   [1] DBI_1.3.0               gson_0.1.0              gridExtra_2.3           rlang_1.1.7            
#   [5] magrittr_2.0.4          otel_0.2.0              compiler_4.5.1          RSQLite_2.4.6          
#   [9] reactome.db_1.92.0      png_0.1-9               systemfonts_1.3.2       vctrs_0.7.2            
#  [13] reshape2_1.4.5          stringr_1.6.0           pkgconfig_2.0.3         crayon_1.5.3           
#  [17] fastmap_1.2.0           XVector_0.50.0          ggraph_2.2.2            UCSC.utils_1.4.0       
#  [21] graph_1.88.1            purrr_1.2.2             bit_4.6.0               xfun_0.55              
#  [25] cachem_1.1.0            graphite_1.54.0         aplot_0.2.9             GenomeInfoDb_1.44.3    
#  [29] jsonlite_2.0.0          blob_1.3.0              tidydr_0.0.6            tweenr_2.0.3           
#  [33] BiocParallel_1.44.0     cluster_2.1.8.2         parallel_4.5.1          R6_2.6.1               
#  [37] stringi_1.8.7           RColorBrewer_1.1-3      cellranger_1.1.0        GOSemSim_2.36.0        
#  [41] Rcpp_1.1.1              Seqinfo_1.0.0           knitr_1.51              ggtangle_0.1.1         
#  [45] R.utils_2.13.0          Matrix_1.7-3            splines_4.5.1           igraph_2.2.3           
#  [49] tidyselect_1.2.1        viridis_0.6.5           qvalue_2.42.0           rstudioapi_0.18.0      
#  [53] dichromat_2.0-0.1       codetools_0.2-20        lattice_0.22-7          tibble_3.3.1           
#  [57] plyr_1.8.9              withr_3.0.2             treeio_1.34.0           KEGGREST_1.50.0        
#  [61] S7_0.2.1                evaluate_1.0.5          gridGraphics_0.5-1      polyclip_1.10-7        
#  [65] scatterpie_0.2.6        Biostrings_2.78.0       pillar_1.11.1           ggtree_4.0.4           
#  [69] ggfun_0.2.0             scales_1.4.0            tidytree_0.4.7          glue_1.8.0             
#  [73] gdtools_0.5.0           lazyeval_0.2.3          tools_4.5.1             ggnewscale_0.5.2       
#  [77] data.table_1.18.2.1     fgsea_1.36.2            ggiraph_0.9.6           graphlayouts_1.2.3     
#  [81] fs_2.0.1                tidygraph_1.3.1         fastmatch_1.1-8         cowplot_1.2.0          
#  [85] grid_4.5.1              tidyr_1.3.2             ape_5.8-1               GenomeInfoDbData_1.2.14
#  [89] nlme_3.1-168            patchwork_1.3.2         ggforce_0.5.0           cli_3.6.6              
#  [93] rappdirs_0.3.4          fontBitstreamVera_0.1.1 viridisLite_0.4.3       gtable_0.3.6           
#  [97] R.methodsS3_1.8.2       yulab.utils_0.2.4       digest_0.6.39           fontquiver_0.2.1       
# [101] ggrepel_0.9.8           ggplotify_0.1.3         htmlwidgets_1.6.4       farver_2.1.2           
# [105] memoise_2.0.1           htmltools_0.5.9         R.oo_1.27.1             lifecycle_1.0.5        
# [109] httr_1.4.8              GO.db_3.22.0            fontLiberation_0.1.0    bit64_4.6.0-1          
# [113] MASS_7.3-65            
