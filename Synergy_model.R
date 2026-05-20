#' ---
#' title: "Synergy DESEQ model"
#' output: html_document
#' date: "2025-11-25"
#' ---
#' 
## ----setup, include=FALSE-------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

#' 
## ----Libraries------------------------------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(DESeq2);library(dplyr);library(tibble);library(readxl);library(writexl)
})

#' 
#' ##------------------------------------------
#' ## DESEQ2 interaction model to test whether the Combo treatment has an effect beyond the individual effects of BMS and Ribo.
#' 
#' # Ribo = R
#' # BMS = B
#' # DMSO = D
#' # Combo = C
#' # Effect of BMS is (B – D)
#' # Effect of Ribo is (R – D)
#' # Effect of Combo is (C – D)
#' 
#' # If C =~ B + R
#'           C – D =~ (B – D) + (R – D)
#'           C – D – (B –D) – (R – D) = C – B – R + D =~ 0
#' # So, If Log2FC > 0 = Synergistically upregulated
#'         Log2FC < 0 =  Synergistically downregulated
#' #-------------------------------------------
#' 
#' 
## ----loading files and preprocess-----------------------------------------------------------------------------------
mat<- read_xlsx("edgeR_new_exploration/deg/gene_filt_counts.xlsx") %>% column_to_rownames("Gene")
metadata<- read.csv("edgeR_new_exploration/metadata.csv")

## Creating metadata variables
metadata <- metadata %>% mutate(
    BMS_treat  = ifelse(Condition %in% c("BMS", "Combo"), 1, 0),
    Ribo_treat = ifelse(Condition %in% c("Ribo", "Combo"), 1, 0))
metadata$BMS_treat  <- factor(metadata$BMS_treat)
metadata$Ribo_treat <- factor(metadata$Ribo_treat)
table(metadata$Condition, metadata$BMS_treat, metadata$Ribo_treat)

#' 
## ----DESEQ model----------------------------------------------------------------------------------------------------
# Creating DESEq2 object 
dds2 <- DESeqDataSetFromMatrix(countData = mat,
                             colData = metadata,
                             design = ~ BMS_treat * Ribo_treat) # ~ BMS_treat * Ribo_treat 
counts(dds2)

## Estimating factors
dds2<-estimateSizeFactors(dds2)
## Performing DESEq
dds2<-DESeq(dds2, full = design(dds2),betaPrior = FALSE)
resultsNames(dds2)
#res<-results(dds, contrast = c("Condition","DMSO","BMS", "Ribo", "Combo"))
res2 = results(dds2, name = "BMS_treat1.Ribo_treat1")
head(res2[order(res2$padj),]) ## order to FDR

### Subsetting the number of DEGs based on 1
#deg<-subset(res, padj < 0.05 & abs(log2FoldChange) > 1 )
CTX_FullTab2 <-as.data.frame(results(dds2, name = "BMS_treat1.Ribo_treat1", cooksCutoff = F, 
                                     independentFiltering = F)) %>% rownames_to_column("Gene")
CTX_DGE2 <- CTX_FullTab2 %>%mutate(Abs = abs(log2FoldChange)) %>%filter(padj < 0.1 & Abs > 0) %>% arrange(desc(Abs))

write_xlsx(CTX_FullTab2, "edgeR_new_exploration/deg/Synergy_modeling/CTX_FullTab.xlsx")
write_xlsx(CTX_DGE2, "edgeR_new_exploration/deg/Synergy_modeling/CTX_DGE.xlsx")
save(dds2, file = "edgeR_new_exploration/deg/Synergy_modeling/dds2.RData")

#' 
## ----r session-info-------------------------------------------------------------------------------------------------
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
#  [1] writexl_1.5.4               readxl_1.4.5                tibble_3.3.1                dplyr_1.2.1                
#  [5] DESeq2_1.48.2               SummarizedExperiment_1.38.1 Biobase_2.68.0              MatrixGenerics_1.20.0      
#  [9] matrixStats_1.5.0           GenomicRanges_1.60.0        GenomeInfoDb_1.44.3         IRanges_2.44.0             
# [13] S4Vectors_0.48.0            BiocGenerics_0.54.1         generics_0.1.4             
# 
# loaded via a namespace (and not attached):
#  [1] SparseArray_1.10.8      lattice_0.22-7          magrittr_2.0.4          evaluate_1.0.5         
#  [5] grid_4.5.1              RColorBrewer_1.1-3      cellranger_1.1.0        jsonlite_2.0.0         
#  [9] Matrix_1.7-3            sessioninfo_1.2.3       httr_1.4.8              UCSC.utils_1.4.0       
# [13] scales_1.4.0            codetools_0.2-20        abind_1.4-8             cli_3.6.6              
# [17] rlang_1.1.7             XVector_0.50.0          DelayedArray_0.36.0     otel_0.2.0             
# [21] S4Arrays_1.10.1         tools_4.5.1             parallel_4.5.1          BiocParallel_1.44.0    
# [25] ggplot2_4.0.2           locfit_1.5-9.12         GenomeInfoDbData_1.2.14 vctrs_0.7.2            
# [29] R6_2.6.1                lifecycle_1.0.5         pkgconfig_2.0.3         pillar_1.11.1          
# [33] gtable_0.3.6            glue_1.8.0              Rcpp_1.1.1              xfun_0.55              
# [37] tidyselect_1.2.1        rstudioapi_0.18.0       knitr_1.51              dichromat_2.0-0.1      
# [41] farver_2.1.2            compiler_4.5.1          S7_0.2.1               
