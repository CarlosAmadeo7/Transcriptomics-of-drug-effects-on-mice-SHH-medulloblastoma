#' ---
#' title: "Pathway activity and Transcription factor inference"
#' output: html_document
#' date: "2025-11-25"
#' ---
#' 
## ----setup, include=FALSE------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

#' 
## ----Libraries-----------------------------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(DESeq2);library(decoupleR);library(dplyr);library(tidyr);library(tibble);library(ggplot2)
  ;library(pheatmap);library(viridis);library(RColorBrewer);library(ggrepel);library(writexl);library(readxl)
  ;library(sessioninfo)
})

#' 
#' #---------------------------------------------
#' # Loading dds file using the DMSO reference
## ------------------------------------------------------------------------------------------------------------------
load("edgeR_new_exploration/deg/dds.RData")

# Using the normalized vst scores from DESEQ2 model
vsd <- vst(dds, blind = FALSE)
vst_mat <- assay(vsd)

### Downloading the get progeny from mice
net <- decoupleR::get_progeny(organism = 'mouse', top = 500)
#net1 <- progeny::getModel(organism = "Mouse", top = 500, decoupleR = TRUE)
net
write_xlsx(net, "edgeR_new_exploration/deg/decoupleR/get_progeny_mice.xlsx")

## Activity inference and progeny model
## Running mlm
sample_acts <- decoupleR::run_mlm(mat = vst_mat, net = net, .source = 'source',.target = 'target',.mor = 'weight', minsize = 5)
sample_acts_mat <- sample_acts %>%
                   tidyr::pivot_wider(id_cols = 'condition', names_from = 'source',values_from = 'score') %>%
                   tibble::column_to_rownames('condition') %>% as.matrix()

# Scale per feature
sample_acts_mat <- scale(sample_acts_mat)
## viz
pal <- viridis::magma(200)
m <- max(abs(sample_acts_mat), na.rm = TRUE)
breaks <- seq(-m, m, length.out = length(pal) + 1)
pdf("edgeR_new_exploration/deg/decoupleR/decoupleR_Progeny_per_sample.pdf", height = 5, width = 10)
pheatmap::pheatmap(mat = t(sample_acts_mat),   fontsize = 14,fontsize_row = 12,fontsize_col = 12,
  color = pal,breaks = breaks,angle_col = 45,treeheight_col = 0,treeheight_row = 8,
  border_color = "grey60",cellwidth = 45,
  cellheight = 15,display_numbers = TRUE,fontsize_number = 13,number_color = "white")
dev.off()

## PROGENy means/ condition
#------------------------------
sample_acts_mat
df_means <- as.data.frame(sample_acts_mat) %>%
  rownames_to_column("sample") %>%
  mutate(Condition = sub("_.*", "", sample),Condition = factor(Condition, levels = c("DMSO", "Ribo", "BMS", "Combo")))
mean_per_condition <- df_means %>%dplyr::select(-sample) %>%
  group_by(Condition) %>%summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>%
  tibble::column_to_rownames("Condition")
## viz
pal <- viridis::magma(200)
m <- max(abs(mean_per_condition), na.rm = TRUE)
breaks <- seq(-m, m, length.out = length(pal) + 1)
pdf("edgeR_new_exploration/deg/decoupleR/decoupleR_Progeny_mean.pdf", height = 4, width = 6)
pheatmap::pheatmap(
  mat = t(mean_per_condition),cluster_cols = F,fontsize = 14,
  fontsize_row = 12,fontsize_col = 12,color = pal,breaks = breaks,angle_col = 45,
  treeheight_col = 0,treeheight_row = 8,border_color = "grey60",cellwidth = 45,
  cellheight = 15,display_numbers = TRUE,fontsize_number = 13,number_color = "white")
dev.off()

## Inference of pathway activities from UNIQUE stat values in Combo vs DMSO
#------------------------------------------------------------------
resultsNames(dds)
res <- results(dds, name = "Condition_Combo_vs_DMSO")

## loading unique genes 
unique <- read_xlsx("edgeR_new_exploration/deg/Unique_combo/Combo_unique_CTX_DGE_shrink.xlsx")
unique_genes <- unique(unique$Gene)
unique_genes <- intersect(unique_genes, rownames(res))
res_unique <- res[unique_genes, ]
res_unique <- res_unique[order(res_unique$padj), ]

## Adjust matrix
deg <- as.data.frame(res_unique) %>%rownames_to_column("ID") %>%
  dplyr::select(ID, stat) %>%filter(!is.na(stat)) %>%column_to_rownames("ID") %>%as.matrix()

contrast_acts <- decoupleR::run_mlm(mat = deg, net = net, .source = 'source', .target = 'target',
                                    .mor = 'weight', minsize = 5)
## Viz
colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")[c(2, 10)])
p <- ggplot2::ggplot(data = contrast_acts, mapping = ggplot2::aes(x = stats::reorder(source, score), y = score)) + 
     ggplot2::geom_bar(mapping = ggplot2::aes(fill = score),color = "black",stat = "identity") +
     ggplot2::scale_fill_gradient2(low = colors[1], mid = "whitesmoke", high = colors[2], midpoint = 0) + 
     ggplot2::theme_minimal() +
     ggplot2::theme(axis.title = element_text(face = "bold", size = 12),
              axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
              axis.text.y = ggplot2::element_text(size = 10, face = "bold"),
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank()) +
     ggplot2::xlab("Pathways")

pdf("edgeR_new_exploration/deg/decoupleR/Progeny/CombovsDMSO/unique_decoupleR_barplot.pdf", height =5, width = 7)
p
dev.off()

write_xlsx(contrast_acts, "edgeR_new_exploration/deg/decoupleR/Progeny/CombovsDMSO/UniquedecoupleR_progeny_CombovsDMSO.xlsx")


## Scatter plot visualization....
pathway <- c('TNFa')
df <- net %>%dplyr::filter(source == pathway) %>%
      dplyr::arrange(target) %>%dplyr::mutate(ID = target, color = "3") %>%tibble::column_to_rownames('target')
inter <- sort(dplyr::intersect(rownames(deg), rownames(df)))
df <- df[inter, ]
df['stat'] <- deg[inter, ]
df <- df %>% dplyr::mutate(color = dplyr::if_else(weight > 0 & stat > 0, '1', color)) %>%
      dplyr::mutate(color = dplyr::if_else(weight > 0 & stat < 0, '2', color)) %>%
      dplyr::mutate(color = dplyr::if_else(weight < 0 & stat > 0, '2', color)) %>%
      dplyr::mutate(color = dplyr::if_else(weight < 0 & stat < 0, '1', color))
colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")[c(2, 10)])
p1 <- ggplot2::ggplot(data = df, mapping = ggplot2::aes(x = weight, y = stat, color = color)) + 
     ggplot2::geom_point(size = 2.5, color = "black") + 
     ggplot2::geom_point(size = 1.5) +
     ggplot2::scale_colour_manual(values = c(colors[2], colors[1], "grey")) +
     ggrepel::geom_label_repel(mapping = ggplot2::aes(label = ID)) + 
     ggplot2::theme_minimal() +
     ggplot2::theme(legend.position = "none") +
     ggplot2::geom_vline(xintercept = 0, linetype = 'dotted') +
     ggplot2::geom_hline(yintercept = 0, linetype = 'dotted') +
     ggplot2::ggtitle(pathway)

pdf("edgeR_new_exploration/deg/decoupleR/Progeny/CombovsDMSO/unique_decoupleR_TGFb.pdf", height =5, width = 7)
p1
dev.off()

#' 
#' # decoupleR for Transcription factors inference
#' #---------------------------------------------
## ------------------------------------------------------------------------------------------------------------------
## Transcription factors 
net <- decoupleR::get_collectri(organism = 'mouse', split_complexes = FALSE)
# Run ulm
sample_acts <- decoupleR::run_ulm(mat = vst_mat, net = net, .source = 'source', 
                                  .target = 'target',.mor = 'mor',  minsize = 5)

n_tfs <- 25 # number of TFs to show
# Transform to wide matrix
sample_acts_mat <- sample_acts %>%tidyr::pivot_wider(id_cols = 'condition', names_from = 'source',values_from = 'score') %>%tibble::column_to_rownames('condition') %>%as.matrix()

# Get top tfs with more variable means across clusters
tfs <- sample_acts %>%dplyr::group_by(source) %>%dplyr::summarise(std = stats::sd(score)) %>%
       dplyr::arrange(-abs(std)) %>%head(n_tfs) %>%dplyr::pull(source)
sample_acts_mat <- sample_acts_mat[,tfs]

# Scale per sample
sample_acts_mat <- scale(sample_acts_mat)
# Choose color palette
colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))
colors.use <- grDevices::colorRampPalette(colors = colors)(100)
my_breaks <- c(seq(-2, 0, length.out = ceiling(100 / 2) + 1),
               seq(0.05, 2, length.out = floor(100 / 2)))

# Plot
pheatmap::pheatmap(mat = sample_acts_mat,color = colors.use,border_color = "white",
                   breaks = my_breaks,cellwidth = 15,cellheight = 15,treeheight_row = 20,treeheight_col = 20)

##-----
# Running contrasts in unique Combo vs DMSO genes
resultsNames(dds)
res <- results(dds, name = "Condition_Combo_vs_DMSO")
## Filtering res based to just get unique ones
unique_genes <- unique(unique$Gene)
unique_genes <- intersect(unique_genes, rownames(res))
res_unique <- res[unique_genes, ]
res_unique <- res_unique[order(res_unique$padj), ]

## Adjust matrix
deg <- as.data.frame(res_unique) %>%rownames_to_column("ID") %>%dplyr::select(ID, stat, log2FoldChange, pvalue ) %>%
  filter(!is.na(stat)) %>%column_to_rownames("ID") %>%as.matrix()
contrast_acts <- decoupleR::run_ulm(mat = deg[, 'stat', drop = FALSE], net = net, .source = 'source', 
                                    .target = 'target',.mor='mor', minsize = 5)
#contrast_acts<- contrast_acts |> filter( p_value < 0.5)
## viz 
# Filter top TFs in both signs
f_contrast_acts <- contrast_acts %>%dplyr::mutate(rnk = NA)
msk <- f_contrast_acts$score > 0
f_contrast_acts[msk, 'rnk'] <- rank(-f_contrast_acts[msk, 'score'])
f_contrast_acts[!msk, 'rnk'] <- rank(-abs(f_contrast_acts[!msk, 'score']))

tfs <- f_contrast_acts %>%dplyr::arrange(rnk) %>%head(n_tfs) %>%dplyr::pull(source)
f_contrast_acts <- f_contrast_acts %>%filter(source %in% tfs)

colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")[c(2, 10)])
p <- ggplot2::ggplot(data = f_contrast_acts, mapping = ggplot2::aes(x = stats::reorder(source, score), y = score)) + 
     ggplot2::geom_bar(mapping = ggplot2::aes(fill = score),color = "black",stat = "identity") +
     ggplot2::scale_fill_gradient2(low = colors[1], mid = "whitesmoke", high = colors[2],  midpoint = 0) + 
     ggplot2::theme_minimal() +
     ggplot2::theme(axis.title = element_text(face = "bold", size = 12),
                    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
              axis.text.y = ggplot2::element_text(size = 10, face = "bold"),
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank()) +
     ggplot2::xlab("TFs")

## Gene-TF viz
tf <- 'Nr4a1' #<- exmaple of TF

df <- net %>%dplyr::filter(source == tf) %>%dplyr::arrange(target) %>%dplyr::mutate(ID = target, color = "3") %>%tibble::column_to_rownames('target')

inter <- sort(dplyr::intersect(rownames(deg), rownames(df)))
df <- df[inter, ]
df$logfc   <- deg[inter, "log2FoldChange"]
df$stat    <- deg[inter, "stat"]
df$p_value <- deg[inter, "pvalue"]

df <- df %>%dplyr::mutate(color = dplyr::if_else(mor > 0 & stat > 0, '1', color)) %>%
      dplyr::mutate(color = dplyr::if_else(mor > 0 & stat < 0, '2', color)) %>%
      dplyr::mutate(color = dplyr::if_else(mor < 0 & stat > 0, '2', color)) %>%
      dplyr::mutate(color = dplyr::if_else(mor < 0 & stat < 0, '1', color))

colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")[c(2, 10)])
p <- ggplot2::ggplot(data = df, mapping = ggplot2::aes(x = logfc, y = -log10(p_value), color = color,size = abs(mor))) + 
     ggplot2::geom_point(size = 2.5, color = "black") + 
     ggplot2::geom_point(size = 1.5) +
     ggplot2::scale_colour_manual(values = c(colors[2], colors[1], "grey")) +
     ggrepel::geom_label_repel(mapping = ggplot2::aes(label = ID,size = 1)) + 
     ggplot2::theme_minimal() +
     ggplot2::theme(legend.position = "none") +
     ggplot2::geom_vline(xintercept = 0, linetype = 'dotted') +
     ggplot2::geom_hline(yintercept = 0, linetype = 'dotted') +
     ggplot2::ggtitle(tf)

#' 
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
#  [1] sessioninfo_1.2.3           readxl_1.4.5                writexl_1.5.4               ggrepel_0.9.8              
#  [5] RColorBrewer_1.1-3          viridis_0.6.5               viridisLite_0.4.3           pheatmap_1.0.13            
#  [9] ggplot2_4.0.2               tibble_3.3.1                tidyr_1.3.2                 dplyr_1.2.1                
# [13] decoupleR_2.14.0            DESeq2_1.48.2               SummarizedExperiment_1.38.1 Biobase_2.68.0             
# [17] MatrixGenerics_1.20.0       matrixStats_1.5.0           GenomicRanges_1.60.0        GenomeInfoDb_1.44.3        
# [21] IRanges_2.44.0              S4Vectors_0.48.0            BiocGenerics_0.54.1         generics_0.1.4             
# 
# loaded via a namespace (and not attached):
#  [1] gtable_0.3.6            xfun_0.55               lattice_0.22-7          vctrs_0.7.2            
#  [5] tools_4.5.1             parallel_4.5.1          pkgconfig_2.0.3         Matrix_1.7-3           
#  [9] S7_0.2.1                lifecycle_1.0.5         GenomeInfoDbData_1.2.14 compiler_4.5.1         
# [13] farver_2.1.2            stringr_1.6.0           codetools_0.2-20        pillar_1.11.1          
# [17] BiocParallel_1.44.0     DelayedArray_0.36.0     abind_1.4-8             parallelly_1.46.1      
# [21] tidyselect_1.2.1        locfit_1.5-9.12         stringi_1.8.7           purrr_1.2.2            
# [25] grid_4.5.1              cli_3.6.6               SparseArray_1.10.8      magrittr_2.0.4         
# [29] S4Arrays_1.10.1         dichromat_2.0-0.1       withr_3.0.2             scales_1.4.0           
# [33] UCSC.utils_1.4.0        XVector_0.50.0          httr_1.4.8              otel_0.2.0             
# [37] gridExtra_2.3           cellranger_1.1.0        evaluate_1.0.5          knitr_1.51             
# [41] rlang_1.1.7             Rcpp_1.1.1              glue_1.8.0              rstudioapi_0.18.0      
# [45] jsonlite_2.0.0          R6_2.6.1               
