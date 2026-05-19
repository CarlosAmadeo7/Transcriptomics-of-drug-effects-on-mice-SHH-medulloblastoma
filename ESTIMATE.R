#' ---
#' title: "ESTIMATE algorithm"
#' output: html_document
#' date: "2025-11-25"
#' ---
#' 
## ----setup, include=FALSE-------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
#' 
## ----Libraries------------------------------------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(readxl);library(dplyr);library(tidyr);library(tibble); library(edgeR); 
  library(immunedeconv);library(ggplot2) ;library(ggpubr); library(Seurat); library(writexl); library(sessioninfo)
})
#' 
## ----loading raw counts---------------------------------------------------------------------------------------------------
counts <- read_xlsx("gene_filt_counts.xlsx")
counts <- counts %>% tibble::column_to_rownames("Gene") %>% as.matrix()
### Converting for ESTIMATE processing 
prepare_expression_data <- function(count_matrix) {
  # 1. Convert to matrix and clean
  if (!is.matrix(count_matrix)) {
    count_matrix <- as.matrix(count_matrix)
  }
  # 2. Remove NA and non-numeric values
  count_matrix[is.na(count_matrix)] <- 0
  count_matrix[!is.numeric(count_matrix)] <- 0
  # 3. Handle duplicate genes (take mean)
  if (any(duplicated(rownames(count_matrix)))) {
    count_matrix <- aggregate(count_matrix, 
                             by = list(rownames(count_matrix)), 
                             FUN = mean)
    rownames(count_matrix) <- count_matrix$Group.1
    count_matrix$Group.1 <- NULL
    count_matrix <- as.matrix(count_matrix)
  }
  # 4. Convert to TPM (required for most methods)
  dge <- DGEList(counts = count_matrix)
  dge <- calcNormFactors(dge)
  tpm_matrix <- cpm(dge, log = FALSE)
  # 5. Ensure proper gene symbols
  rownames(tpm_matrix) <- gsub("\\..*", "", rownames(tpm_matrix))
  return(tpm_matrix)
}
matrix <- prepare_expression_data(counts)
#' 
## ----Converting to human and runnig ESTIMATE------------------------------------------------------------------------------
dataset_humanGenes <- convert_human_mouse_genes(matrix, convert_to = 'human')
# Running ESTIMATE
ESTIMATE<- immunedeconv::deconvolute_estimate(dataset_humanGenes)
# Viz
ESTIMATE<- ESTIMATE %>% tibble::rownames_to_column("cell_type")
ESTIMATE_long <- ESTIMATE %>%pivot_longer(cols = -cell_type,names_to = "sample",values_to = "score") %>% 
  mutate(condition = sub("_.*", "", sample))
ESTIMATE_long$condition <- factor(ESTIMATE_long$condition,
  levels = c("DMSO", "BMS", "Ribo", "Combo"),
  labels = c("Vehicle", "BMS-986158", "Ribociclib", "BMS-986158/Ribociclib")
)
## Example colors: #4DBBD5","#E64B35", "#00A087", "#3C5488"
pdf("deconvolution/ESTIMATES_plot.pdf", height = 5, width = 12)
ggpubr::ggboxplot(
  ESTIMATE_long,x = "condition",y = "score",
  color = "condition", add = "jitter", legend = "right", palette = c(
    "Vehicle" = "#D8CE7A","BMS-986158"= "#A5D996","Ribociclib"  = "#8259D9","BMS-986158/Ribociclib" = "blue4")) +
  facet_wrap(~cell_type,
             scales = "free_y",
             ncol = 4) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    axis.title.x = element_blank(),
    strip.text = element_text(size = 14, face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 11)
  ) + RotatedAxis()
dev.off()

# Writing ESTIMATE results 
write_xlsx(ESTIMATE, "ESTIMATE_results.xlsx")
#' 
## ----r session-info-------------------------------------------------------------------------------------------------------
sessioninfo::session_info()
#' 
#' 

# R version 4.5.1 (2025-06-13 ucrt)
# Platform: x86_64-w64-mingw32/x64
# Running under: Windows 11 x64 (build 26200)
# 
# Matrix products: default
#   LAPACK version 3.12.1
# 
# locale:
# [1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8    LC_MONETARY=English_United States.utf8
# [4] LC_NUMERIC=C                           LC_TIME=English_United States.utf8    
# 
# time zone: America/New_York
# tzcode source: internal
# 
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#  [1] sessioninfo_1.2.3  writexl_1.5.4      Seurat_5.4.0       SeuratObject_5.3.0 sp_2.2-1           ggpubr_0.6.3      
#  [7] ggplot2_4.0.2      immunedeconv_2.1.0 EPIC_1.1.7         edgeR_4.8.2        limma_3.66.0       tibble_3.3.1      
# [13] tidyr_1.3.2        dplyr_1.2.1        readxl_1.4.5      
# 
# loaded via a namespace (and not attached):
#   [1] RcppAnnoy_0.0.23        splines_4.5.1           later_1.4.8             filelock_1.0.3         
#   [5] cellranger_1.1.0        polyclip_1.10-7         preprocessCore_1.71.2   XML_3.99-0.23          
#   [9] fastDummies_1.7.5       lifecycle_1.0.5         httr2_1.2.2             rstatix_0.7.3          
#  [13] globals_0.19.1          lattice_0.22-7          MASS_7.3-65             backports_1.5.1        
#  [17] magrittr_2.0.4          rmarkdown_2.31          plotly_4.12.0           yaml_2.3.12            
#  [21] httpuv_1.6.17           otel_0.2.0              sctransform_0.4.3       spam_2.11-3            
#  [25] spatstat.sparse_3.1-0   reticulate_1.46.0       cowplot_1.2.0           pbapply_1.7-4          
#  [29] DBI_1.3.0               RColorBrewer_1.1-3      abind_1.4-8             Rtsne_0.17             
#  [33] purrr_1.2.2             BiocGenerics_0.54.1     rappdirs_0.3.4          sva_3.56.0             
#  [37] GenomeInfoDbData_1.2.14 IRanges_2.44.0          S4Vectors_0.48.0        data.tree_1.2.0        
#  [41] ggrepel_0.9.8           irlba_2.3.7             spatstat.utils_3.2-2    listenv_0.10.1         
#  [45] genefilter_1.90.0       goftest_1.2-3           RSpectra_0.16-2         spatstat.random_3.4-5  
#  [49] annotate_1.88.0         fitdistrplus_1.2-6      parallelly_1.46.1       codetools_0.2-20       
#  [53] xml2_1.5.2              tidyselect_1.2.1        UCSC.utils_1.4.0        farver_2.1.2           
#  [57] spatstat.explore_3.8-0  matrixStats_1.5.0       stats4_4.5.1            BiocFileCache_3.0.0    
#  [61] Seqinfo_1.0.0           jsonlite_2.0.0          progressr_0.19.0        Formula_1.2-5          
#  [65] ggridges_0.5.7          survival_3.8-3          tools_4.5.1             progress_1.2.3         
#  [69] ica_1.0-3               Rcpp_1.1.1              glue_1.8.0              gridExtra_2.3          
#  [73] xfun_0.55               mgcv_1.9-3              MatrixGenerics_1.20.0   GenomeInfoDb_1.44.3    
#  [77] withr_3.0.2             ComICS_1.0.4            fastmap_1.2.0           digest_0.6.39          
#  [81] R6_2.6.1                mime_0.13               scattermore_1.2         tensor_1.5.1           
#  [85] spatstat.data_3.1-9     dichromat_2.0-0.1       biomaRt_2.64.0          mMCPcounter_1.1.0      
#  [89] RSQLite_2.4.6           generics_0.1.4          data.table_1.18.2.1     prettyunits_1.2.0      
#  [93] httr_1.4.8              htmlwidgets_1.6.4       uwot_0.2.4              pkgconfig_2.0.3        
#  [97] gtable_0.3.6            blob_1.3.0              lmtest_0.9-40           S7_0.2.1               
# [101] XVector_0.50.0          htmltools_0.5.9         carData_3.0-6           dotCall64_1.2          
# [105] scales_1.4.0            Biobase_2.68.0          png_0.1-9               spatstat.univar_3.1-7  
# [109] knitr_1.51              rstudioapi_0.18.0       reshape2_1.4.5          tzdb_0.5.0             
# [113] nlme_3.1-168            curl_7.0.0              cachem_1.1.0            zoo_1.8-15             
# [117] stringr_1.6.0           testit_0.17             KernSmooth_2.23-26      parallel_4.5.1         
# [121] miniUI_0.1.2            AnnotationDbi_1.72.0    pillar_1.11.1           grid_4.5.1             
# [125] vctrs_0.7.2             RANN_2.6.2              promises_1.5.0          car_3.1-5              
# [129] dbplyr_2.5.2            xtable_1.8-8            cluster_2.1.8.2         evaluate_1.0.5         
# [133] readr_2.2.0             cli_3.6.6               locfit_1.5-9.12         compiler_4.5.1         
# [137] rlang_1.1.7             crayon_1.5.3            future.apply_1.20.2     ggsignif_0.6.4         
# [141] labeling_0.4.3          plyr_1.8.9              stringi_1.8.7           deldir_2.0-4           
# [145] viridisLite_0.4.3       BiocParallel_1.44.0     Biostrings_2.78.0       lazyeval_0.2.3         
# [149] spatstat.geom_3.7-3     Matrix_1.7-3            RcppHNSW_0.6.0          hms_1.1.4              
# [153] patchwork_1.3.2         bit64_4.6.0-1           future_1.70.0           KEGGREST_1.50.0        
# [157] statmod_1.5.1           shiny_1.13.0            ROCR_1.0-12             igraph_2.2.3           
# [161] broom_1.0.12            memoise_2.0.1           bit_4.6.0              
