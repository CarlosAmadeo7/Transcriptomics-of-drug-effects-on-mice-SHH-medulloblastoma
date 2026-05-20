#' ---
#' title: "Analysis on Ribociclib , BMS-986156 , Combo vs Vehicle"
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
  library(writexl);library(scToppR);
  library(clusterProfiler);library(readxl);
  library(ComplexHeatmap); library(ggVennDiagram); library(ggvenn); library(UpSetR);library(sessioninfo)
})

#' 
## ----Loading counts------------------------------------------------------------------------------------------------
gene_counts<- read_xlsx("edgeR_new_exploration/counts.xlsx")
counts<- gene_counts |> filter(gene_biotype == "protein_coding") |>
  dplyr::select(gene_name, 2:13) |> dplyr::select(gene_name,
    DMSO_1, DMSO_2, DMSO_3,
    BMS_2, BMS_3, BMS_4,
    Ribo_1, Ribo_2, Ribo_3,
    Combo_1, Combo_3, Combo_4)
rownames(counts)<- NULL

p_exp <- counts |> dplyr::distinct(gene_name, .keep_all = TRUE)
rownames(p_exp)<- NULL
p_exp <- p_exp |> column_to_rownames('gene_name')

# Reading metadata
metadata<- read.csv("edgeR_new_exploration/metadata.csv")
metadata<- metadata |> column_to_rownames("X")
metadata$Condition<-factor(trimws(metadata$Condition))
metadata$Condition<-factor(metadata$Condition, levels = c("DMSO","BMS", "Ribo", "Combo"))
levels(metadata$Condition)
all(colnames(p_exp) == rownames(metadata))

#' 
## ----CPM normalization---------------------------------------------------------------------------------------------
## Converting to cpms: Counts per million 
cpm <-future_apply(p_exp, 2, function(x) x/sum(as.numeric(x)) * 10^6) 

## Filter low cpm counts for both conditions 
filter<-apply(cpm, 1, function(x) all(x[1:3]>=0.5) | all(x[4:6]>=0.5) | all(x[7:9]>=0.5) | all(x[10:12]>=0.5))
counts_filt <- p_exp[filter,]  ## 13772 genes
cpm_filt <- cpm[filter,]   ## 13772 genes

## Saving cpm and counts filter
cpm_filt_xlsx<-cpm_filt %>% as.data.frame() %>% tibble::rownames_to_column(var = "Gene")
write_xlsx(cpm_filt_xlsx, "edgeR_new_exploration/deg/gene_cpm.xlsx")

counts_filt_xlsx<-counts_filt %>%as.data.frame() %>%tibble::rownames_to_column(var = "Gene")
write_xlsx(counts_filt_xlsx, "edgeR_new_exploration/deg/gene_filt_counts.xlsx")

## Applying log_transformation
logCPM<- log2(cpm_filt +1)
boxplot(logCPM)
## Normalizing quantiles from cpm counts
p<- normalize.quantiles(as.matrix(logCPM))
rownames(p)<-rownames(logCPM)
colnames(p) <-colnames(logCPM)
boxplot(p)
## Performing PCA using the quantile normalization from the CPMs and Visualize
metadata$Subject_ID<- rownames(metadata)
pca.Samples<-prcomp(t(p))
PCi<-data.frame(pca.Samples$x, Condition= metadata$Condition, ID= metadata$Subject_ID )
eig <- (pca.Samples$sdev)^2 ### Calculating the eigenvalue
variance <- eig*100/sum(eig) ## Calculating variance
pdf("edgeR_new_exploration/deg/PCA_cpms.pdf",width=8,height=5,useDingbats=FALSE)
ggscatter(PCi,x = "PC1",y = "PC2",color = "Condition",
  palette = c("#D8CE7A", "#A5D996", "#8259D9", "blue4"),size = 6,ellipse = TRUE,
  ellipse.type = "confidence",ellipse.alpha = 0.2) +
  xlab(paste("PC1 (", round(variance[1],1), "% )")) +
  ylab(paste("PC2 (", round(variance[2],1), "% )")) +
  labs(color = "Condition", fill = "Condition") +
  scale_color_manual(
    values = c("#D8CE7A", "#A5D996", "#8259D9", "blue4"),
    labels = c("Vehicle", "BMS-986158", "Ribociclib", "BMS-986158/Ribociclib")
  ) +
  scale_fill_manual(
    values = c("#D8CE7A", "#A5D996", "#8259D9", "blue4"),
    labels = c("Vehicle", "BMS-986158", "Ribociclib", "BMS-986158/Ribociclib")
  ) +
  theme_classic() +
  theme(axis.title = element_text(size = 16),axis.text = element_text(size = 16),
    legend.title = element_text(size = 14),legend.text = element_text(size = 12),
    axis.ticks = element_line(linewidth = 0.8),
    axis.ticks.length = unit(0.2, "cm"))
dev.off()

#' 
#' #------- DESEQ2
#' ## Differential expression analysis using DMSO as reference
## ----Combo vs DMSO-------------------------------------------------------------------------------------------------
## Setting metadata condition as : DMSO as reference
metadata$Condition<- factor(metadata$Condition)
metadata$Condition <- relevel(metadata$Condition, ref = "DMSO")

## Creating DESEq object 
dds<- DESeqDataSetFromMatrix(countData = p_exp,colData = metadata, design = ~ Condition)
#counts(dds)
## Estimating factors
dds<-estimateSizeFactors(dds)
## Performing DESEq
dds<-DESeq(dds, full = design(dds),betaPrior = FALSE)
resultsNames(dds)

#--------- Obtaining results
#-----------------------------------
#----------------------------------- COMBO VS DMSO
cd <- "Condition_Combo_vs_DMSO"  
res = results(dds, name = cd)
head(res[order(res$padj),]) ## order to FDR
## Creating tmps
CTX_FullTab <-as.data.frame(results(dds, name = cd, cooksCutoff = F, independentFiltering = F))%>% rownames_to_column("Gene")
CTX_DGE <- CTX_FullTab %>%mutate(Abs = abs(log2FoldChange)) %>%filter(padj < 0.05 & Abs > 0.3) %>% arrange(desc(Abs))

## Establishing filepath --< Changing according to purposes
file_path<- "edgeR_new_exploration/Refinment/deg/Combo_DMSO/"
## Saving
write_xlsx(CTX_FullTab, file.path(file_path, "CTX_FullTab.xlsx"))
write_xlsx(CTX_DGE, file.path(file_path, "CTX_DGE.xlsx"))

#--
## Shrinkage
res_shrink<-lfcShrink(dds, coef = cd, type = "ashr" )
degs_shrink<- subset(res_shrink, padj < 0.05 & abs(log2FoldChange) > 0.3 )

## Creating tmps
CTX_FullTab_shrink<- as.data.frame(res_shrink) %>%rownames_to_column("Gene")
CTX_DGE_shrink<-as.data.frame(degs_shrink) %>%rownames_to_column("Gene")
## Saving
write_xlsx(CTX_FullTab_shrink, file.path(file_path, "CTX_FullTab_shrink.xlsx"))
write_xlsx(CTX_DGE_shrink, file.path(file_path, "CTX_DGE_shrink.xlsx"))

save(CTX_FullTab, CTX_DGE, CTX_FullTab_shrink, CTX_DGE_shrink, metadata, metadata, 
     file = file.path(file_path, "CTX_Dge_Data.RData"))
save(dds, file = file.path(file_path, "Ref_DMSO.RData"))

#' 
#----------------------------------- RIBO VS DMSO 
cd <- "Condition_Ribo_vs_DMSO"  
res = results(dds, name = cd)
head(res[order(res$padj),]) ## order to FDR
## Saving
CTX_FullTab <-as.data.frame(results(dds, name = cd, cooksCutoff = F, independentFiltering = F))%>% rownames_to_column("Gene")
CTX_DGE <- CTX_FullTab %>%mutate(Abs = abs(log2FoldChange)) %>%filter(padj < 0.05 & Abs > 0.3) %>% arrange(desc(Abs))

## Establishing filepath --< Changing accoridng to purposes
file_path<- "edgeR_new_exploration/Refinment/deg/Ribo_DMSO/"

write_xlsx(CTX_FullTab, file.path(file_path, "CTX_FullTab.xlsx"))
write_xlsx(CTX_DGE, file.path(file_path, "CTX_DGE.xlsx"))

#--
## Shrinkage
res_shrink<-lfcShrink(dds, coef = cd, type = "ashr" )
degs_shrink<- subset(res_shrink, padj < 0.05 & abs(log2FoldChange) > 0.3 )

### Creating tmps
CTX_FullTab_shrink<- as.data.frame(res_shrink) %>%rownames_to_column("Gene")
CTX_DGE_shrink<-as.data.frame(degs_shrink) %>%rownames_to_column("Gene")

write_xlsx(CTX_FullTab_shrink, file.path(file_path, "CTX_FullTab_shrink.xlsx"))
write_xlsx(CTX_DGE_shrink, file.path(file_path, "CTX_DGE_shrink.xlsx"))
save(CTX_FullTab, CTX_DGE, CTX_FullTab_shrink, CTX_DGE_shrink, metadata, metadata, 
     file = file.path(file_path, "CTX_Dge_Data.RData"))

#----------------------------------- BMS VS DMSO 
cd <- "Condition_BMS_vs_DMSO"  
res = results(dds, name = cd)
head(res[order(res$padj),]) ## order to FDR
#### Saving
CTX_FullTab <-as.data.frame(results(dds, name = cd, cooksCutoff = F, independentFiltering = F))%>% rownames_to_column("Gene")
CTX_DGE <- CTX_FullTab %>%
  mutate(Abs = abs(log2FoldChange)) %>%
  filter(padj < 0.05 & Abs > 0.3) %>% 
  arrange(desc(Abs))
## Establishing filepath --< Changing accoridng to purposes
file_path<- "edgeR_new_exploration/Refinment/deg/BMS_DMSO/"

write_xlsx(CTX_FullTab, file.path(file_path, "CTX_FullTab.xlsx"))
write_xlsx(CTX_DGE, file.path(file_path, "CTX_DGE.xlsx"))

#--
## Shrinkage
res_shrink<-lfcShrink(dds, coef = cd, type = "ashr" )
degs_shrink<- subset(res_shrink, padj < 0.05 & abs(log2FoldChange) > 0.3 )

### Creating tmp1
CTX_FullTab_shrink<- as.data.frame(res_shrink) %>%rownames_to_column("Gene")
CTX_DGE_shrink<-as.data.frame(degs_shrink) %>%rownames_to_column("Gene")

write_xlsx(CTX_FullTab_shrink, file.path(file_path, "CTX_FullTab_shrink.xlsx"))
write_xlsx(CTX_DGE_shrink, file.path(file_path, "CTX_DGE_shrink.xlsx"))
save(CTX_FullTab, CTX_DGE, CTX_FullTab_shrink, CTX_DGE_shrink, metadata, metadata, file = file.path(file_path, "CTX_Dge_Data.RData"))

#' 
#' ## Differential expression analysis using BMS as reference
## ----Combo vs BMS----
## Setting metadata condition as : BMS as reference
metadata$Condition<- factor(metadata$Condition)
metadata$Condition <- relevel(metadata$Condition, ref = "BMS")

## Creating DESEq object 
dds<- DESeqDataSetFromMatrix(countData = p_exp,colData = metadata, design = ~ Condition)
#counts(dds)
## Estimating factors
dds<-estimateSizeFactors(dds)
## Performing DESEq
dds<-DESeq(dds, full = design(dds),betaPrior = FALSE)
resultsNames(dds)

#----------------------------------- COMBO VS BMS
cd <- "Condition_Combo_vs_BMS"
res = results(dds, name = cd)
head(res[order(res$padj),]) ## order to FDR
## Creating tmps
CTX_FullTab <-as.data.frame(results(dds, name = cd, cooksCutoff = F, independentFiltering = F))%>% rownames_to_column("Gene")
CTX_DGE <- CTX_FullTab %>%mutate(Abs = abs(log2FoldChange)) %>%filter(padj < 0.05 & Abs > 0.3) %>% arrange(desc(Abs))

## Establishing filepath --< Changing according to purposes
file_path<- "edgeR_new_exploration/Refinment/deg/Combo_BMS/"
## Saving
write_xlsx(CTX_FullTab, file.path(file_path, "CTX_FullTab.xlsx"))
write_xlsx(CTX_DGE, file.path(file_path, "CTX_DGE.xlsx"))

#--
## Shrinkage
res_shrink<-lfcShrink(dds, coef = cd, type = "ashr" )
degs_shrink<- subset(res_shrink, padj < 0.05 & abs(log2FoldChange) > 0.3 )

## Creating tmps
CTX_FullTab_shrink<- as.data.frame(res_shrink) %>%rownames_to_column("Gene")
CTX_DGE_shrink<-as.data.frame(degs_shrink) %>%rownames_to_column("Gene")
## Saving
write_xlsx(CTX_FullTab_shrink, file.path(file_path, "CTX_FullTab_shrink.xlsx"))
write_xlsx(CTX_DGE_shrink, file.path(file_path, "CTX_DGE_shrink.xlsx"))

save(CTX_FullTab, CTX_DGE, CTX_FullTab_shrink, CTX_DGE_shrink, metadata, metadata, file = file.path(file_path, "CTX_Dge_Data.RData"))
save(dds, file = file.path(file_path, "BMS_dds.RData"))

#' 
#----------------------------------- RIBO VS BMS 
cd <- "Condition_Ribo_vs_BMS"  
res = results(dds, name = cd)
head(res[order(res$padj),]) ## order to FDR
## Saving
CTX_FullTab <-as.data.frame(results(dds, name = cd, cooksCutoff = F, independentFiltering = F))%>% rownames_to_column("Gene")
CTX_DGE <- CTX_FullTab %>%mutate(Abs = abs(log2FoldChange)) %>%filter(padj < 0.05 & Abs > 0.3) %>% arrange(desc(Abs))

## Establishing filepath --< Changing accoridng to purposes
file_path<- "edgeR_new_exploration/Refinment/deg/Ribo_BMS/"

write_xlsx(CTX_FullTab, file.path(file_path, "CTX_FullTab.xlsx"))
write_xlsx(CTX_DGE, file.path(file_path, "CTX_DGE.xlsx"))

#---
## Shrinkage
res_shrink<-lfcShrink(dds, coef = cd, type = "ashr" )
degs_shrink<- subset(res_shrink, padj < 0.05 & abs(log2FoldChange) > 0.3 )

### Creating tmps
CTX_FullTab_shrink<- as.data.frame(res_shrink) %>%rownames_to_column("Gene")
CTX_DGE_shrink<-as.data.frame(degs_shrink) %>%rownames_to_column("Gene")

write_xlsx(CTX_FullTab_shrink, file.path(file_path, "CTX_FullTab_shrink.xlsx"))
write_xlsx(CTX_DGE_shrink, file.path(file_path, "CTX_DGE_shrink.xlsx"))
save(CTX_FullTab, CTX_DGE, CTX_FullTab_shrink, CTX_DGE_shrink, metadata, metadata, file = file.path(file_path, "CTX_Dge_Data.RData"))

#' 
#' ## Differential expression analysis using Ribo as reference
## ----Combo vs Ribo-------------------------------------------------------------------------------------------------
## Setting metadata condition as : Ribo as reference
metadata$Condition<- factor(metadata$Condition)
metadata$Condition <- relevel(metadata$Condition, ref = "Ribo")

## Creating DESEq object 
dds<- DESeqDataSetFromMatrix(countData = p_exp,colData = metadata, design = ~ Condition)
#counts(dds)
## Estimating factors
dds<-estimateSizeFactors(dds)
## Performing DESEq
dds<-DESeq(dds, full = design(dds),betaPrior = FALSE)
resultsNames(dds)

#----------------------------------- COMBO VS RIBO
cd <- "Condition_Combo_vs_Ribo"  
res = results(dds, name = cd)
head(res[order(res$padj),]) ## order to FDR
## Creating tmps
CTX_FullTab <-as.data.frame(results(dds, name = cd, cooksCutoff = F, independentFiltering = F))%>% rownames_to_column("Gene")
CTX_DGE <- CTX_FullTab %>%mutate(Abs = abs(log2FoldChange)) %>%filter(padj < 0.05 & Abs > 0.3) %>% arrange(desc(Abs))

## Establishing filepath --< Changing according to purposes
file_path<- "edgeR_new_exploration/Refinment/deg/Combo_Ribo/"
## Saving
write_xlsx(CTX_FullTab, file.path(file_path, "CTX_FullTab.xlsx"))
write_xlsx(CTX_DGE, file.path(file_path, "CTX_DGE.xlsx"))

#--
## Shrinkage
res_shrink<-lfcShrink(dds, coef = cd, type = "ashr" )
degs_shrink<- subset(res_shrink, padj < 0.05 & abs(log2FoldChange) > 0.3 )

## Creating tmps
CTX_FullTab_shrink<- as.data.frame(res_shrink) %>%rownames_to_column("Gene")
CTX_DGE_shrink<-as.data.frame(degs_shrink) %>%rownames_to_column("Gene")
## Saving
write_xlsx(CTX_FullTab_shrink, file.path(file_path, "CTX_FullTab_shrink.xlsx"))
write_xlsx(CTX_DGE_shrink, file.path(file_path, "CTX_DGE_shrink.xlsx"))

save(CTX_FullTab, CTX_DGE, CTX_FullTab_shrink, CTX_DGE_shrink, metadata, metadata, file = file.path(file_path, "CTX_Dge_Data.RData"))
save(dds, file = file.path(file_path, "Ribo_dds.RData"))


#' 
#' # Computing Unique genes/condition
## ------------------------------------------------------------------------------------------------------------------
# Loading DEGs dataframes
Ribo<- read_xlsx("edgeR_new_exploration/deg/Ribo_DMSO/CTX_DGE_shrink.xlsx")
BMS<- read_xlsx("edgeR_new_exploration/deg/BMS_DMSO/CTX_DGE_shrink.xlsx")
Combo<- read_xlsx("edgeR_new_exploration/deg/Combo_DMSO/CTX_DGE_shrink.xlsx")
Combo_BMS<- read_xlsx("edgeR_new_exploration/deg/Combo_BMS/CTX_DGE_shrink.xlsx")
Combo_Ribo<- read_xlsx("edgeR_new_exploration/deg/Combo_Ribo/CTX_DGE_shrink.xlsx")

# Creating temporary dfs
RiboDEGs <- Ribo %>%
  mutate(LOG = -log10(padj), ABS = abs(log2FoldChange)) %>%
  mutate(Threshold = if_else(padj < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
  mutate(Direction = case_when(log2FoldChange > 0.3 & padj < 0.05 ~ "UpReg", log2FoldChange < -0.3 & padj < 0.05 ~ "DownReg")) 
BMSDEGs <- BMS %>%
  mutate(LOG = -log10(padj), ABS = abs(log2FoldChange)) %>%
  mutate(Threshold = if_else(padj < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
  mutate(Direction = case_when(log2FoldChange > 0.3 & padj < 0.05 ~ "UpReg", log2FoldChange < -0.3 & padj < 0.05 ~ "DownReg")) 
ComboDEGs <- Combo %>%
  mutate(LOG = -log10(padj), ABS = abs(log2FoldChange)) %>%
  mutate(Threshold = if_else(padj < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
  mutate(Direction = case_when(log2FoldChange > 0.3 & padj < 0.05 ~ "UpReg", log2FoldChange < -0.3 & padj < 0.05 ~ "DownReg")) 
Combo_BMSDEGs <- Combo_BMS %>%
  mutate(LOG = -log10(padj), ABS = abs(log2FoldChange)) %>%
  mutate(Threshold = if_else(padj < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
  mutate(Direction = case_when(log2FoldChange > 0.3 & padj < 0.05 ~ "UpReg", log2FoldChange < -0.3 & padj < 0.05 ~ "DownReg")) 
Combo_RiboDEGs <- Combo_Ribo %>%
  mutate(LOG = -log10(padj), ABS = abs(log2FoldChange)) %>%
  mutate(Threshold = if_else(padj < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
  mutate(Direction = case_when(log2FoldChange > 0.3 & padj < 0.05 ~ "UpReg", log2FoldChange < -0.3 & padj < 0.05 ~ "DownReg")) 


## BMS driven
## Genes significantly changed in both BMS and Combo vs DMSO, not changed by Ribo, and Combo is not significantly different from BMS.
BMS_driven <- intersect(BMSDEGs$Gene, ComboDEGs$Gene) # Intersect BMS and Combo
BMS_driven <- setdiff(BMS_driven, RiboDEGs$Gene) # exclude Ribo-driven
BMS_driven <- setdiff(BMS_driven, Combo_BMS$Gene) # Combo not different from BMS
#BMS_driven <- setdiff(BMS_driven, CTX_DGE2$Gene) # keep synergy separate
write_xlsx(BMS_driven, "edgeR_new_exploration/deg/Unique_BMS/unique_BMS_DEGs.xlsx")

## Ribo driven 
## Genes that change with Ribo and Combo, independent of BMS.
Ribo_driven <- intersect(RiboDEGs$Gene, ComboDEGs$Gene) # Intersect Ribo and Combo
Ribo_driven <- setdiff(Ribo_driven, BMSDEGs$Gene) # exclude BMS-driven
Ribo_driven <- setdiff(Ribo_driven, Combo_RiboDEGs$Gene) # Combo not different from Ribo
#Ribo_driven <- setdiff(Ribo_driven, CTX_DGE2$Gene) # keep synergy separate
write_xlsx(Ribo_driven, "edgeR_new_exploration/deg/Unique_Ribo/unique_Ribo_DEGs.xlsx")

### Combo driven
## Changed only in Combo vs DMSO, not in either single-drug (within thresholds)
Combo_driven <- setdiff(ComboDEGs$Gene, union(BMSDEGs$Gene, RiboDEGs$Gene))
#Combo_driven <- setdiff(Combo_driven, CTX_DGE2$Gene)  # synergy kept separate
write_xlsx(Combo_driven, "edgeR_new_exploration/deg/Unique_combo/Combo_unique_CTX_DGE_shrink.xlsx")


### Unique genes per condition
gene_lists <- list(`BMS-986158` = BMS_driven, Ribociclib= Ribo_driven, `BMS-986158/Ribociclib`= Combo_driven)
p4<-ggvenn(gene_lists, fill_color = c("#A5D996", "#8259D9", "blue4"),
  stroke_size = 0.8, set_name_size = 5, show_percentage = F, text_size = 6,auto_scale = F)
pdf("edgeR_new_exploration/deg/Unique_combo/VennDiagram.pdf", height = 5)
p4
dev.off()

#' 
## ----Overlap VennDiagram-------------------------------------------------------------------------------------------
annotate_degs <- function(df) {df %>% mutate(LOG = -log10(padj),ABS = abs(log2FoldChange),Threshold = padj < 0.05 & ABS > 0.3,
      Direction = case_when(log2FoldChange > 0.3 & padj < 0.05 ~ "UpReg",log2FoldChange < -0.3 & padj < 0.05 ~ "DownReg",TRUE ~ "NotSig"))}

RiboDEGs <- annotate_degs(Ribo)
BMSDEGs <- annotate_degs(BMS)
ComboDEGs <- annotate_degs(Combo)
Combo_BMSDEGs <- annotate_degs(Combo_BMS)
Combo_RiboDEGs <- annotate_degs(Combo_Ribo)

Ribo_sig <- RiboDEGs %>%filter(Threshold) %>% pull(Gene) %>% unique()
BMS_sig <- BMSDEGs %>% filter(Threshold) %>% pull(Gene) %>% unique()
Combo_sig <- ComboDEGs %>% filter(Threshold) %>% pull(Gene) %>% unique()
Combo_BMS_sig <- Combo_BMSDEGs %>% filter(Threshold) %>% pull(Gene) %>% unique()
Combo_Ribo_sig <- Combo_RiboDEGs %>% filter(Threshold) %>% pull(Gene) %>% unique()

gene_lists_overlap <- list(`BMS-986158` = BMS_sig,`Ribociclib` = Ribo_sig,`BMS-986158/Ribociclib` = Combo_sig)
p_overlap <- ggvenn(gene_lists_overlap,fill_color = c("#A5D996", "#8259D9", "blue4"),stroke_size = 0.8,
  set_name_size = 5,show_percentage = FALSE,text_size = 6,auto_scale = FALSE)

pdf("edgeR_new_exploration/deg/Unique_combo/VennDiagram_overlap.pdf", height = 5)
p_overlap
dev.off()

#' 
## ----UpSet plot----------------------------------------------------------------------------------------------------
# Annotate DEG tables
RiboDEGs <- annotate_degs(Ribo)
BMSDEGs <- annotate_degs(BMS)
ComboDEGs <- annotate_degs(Combo)
# Extract significant genes
Ribo_sig <- RiboDEGs %>%filter(Threshold) %>%pull(Gene) %>%unique()
BMS_sig <- BMSDEGs %>%filter(Threshold) %>%pull(Gene) %>%unique()
Combo_sig <- ComboDEGs %>%filter(Threshold) %>%pull(Gene) %>%unique()
# Build gene sets
sets_deg <- list(`BMS-986158` = BMS_sig,`Ribociclib` = Ribo_sig,`BMS-986158/Ribociclib` = Combo_sig)
# Convert to UpSetR format
tmp_upset <- bind_rows(lapply(names(sets_deg), function(nm) {
  data.frame(Gene = sets_deg[[nm]],Class = nm)}))
l <- split(as.character(tmp_upset$Gene), tmp_upset$Class)
# Set desired order
l <- l[c("Ribociclib","BMS-986158","BMS-986158/Ribociclib")]
# Metadata for set colors
metadata <- data.frame(sets = names(l),ToTGene = as.integer(lengths(l)))
# Colors
cols <- c(`BMS-986158` = "#A5D996",Ribociclib = "#8259D9",`BMS-986158/Ribociclib` = "blue4")
# UpSet plot
p <- UpSetR::upset(
  fromList(l),nsets = 3,sets = names(l),keep.order = TRUE,order.by = "freq",
  # Increase text sizes
  text.scale = c(
    1.8, # intersection size title
    1.8, # intersection size numbers
    1.8, # set size title
    1.8, # set size numbers
    1.8, # set names
    1.8  # numbers above bars
  ),
  number.angles = 0,
  set.metadata = list(
    data = metadata,
    plots = list(
      list(type = "matrix_rows",column = "sets",colors = cols,alpha = 0.5))))
pdf("edgeR_new_exploration/deg/UpSetPlot.pdf", height = 5, width = 9)
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
#  [1] sessioninfo_1.2.3           UpSetR_1.4.0                ggvenn_0.1.19               ggVennDiagram_1.5.7        
#  [5] ComplexHeatmap_2.26.0       readxl_1.4.5                clusterProfiler_4.18.4      scToppR_0.99.4             
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
#   [1] fs_2.0.1                enrichplot_1.30.4       httr_1.4.8              doParallel_1.0.17      
#   [5] tools_4.5.1             backports_1.5.1         R6_2.6.1                lazyeval_0.2.3         
#   [9] GetoptLong_1.1.1        withr_3.0.2             gridExtra_2.3           cli_3.6.6              
#  [13] scatterpie_0.2.6        mvtnorm_1.3-6           S7_0.2.1                systemfonts_1.3.2      
#  [17] yulab.utils_0.2.4       gson_0.1.0              DOSE_4.4.0              R.utils_2.13.0         
#  [21] dichromat_2.0-0.1       parallelly_1.46.1       rstudioapi_0.18.0       RSQLite_2.4.6          
#  [25] gridGraphics_0.5-1      shape_1.4.6.1           car_3.1-5               zip_2.3.3              
#  [29] GO.db_3.22.0            Matrix_1.7-3            abind_1.4-8             R.methodsS3_1.8.2      
#  [33] lifecycle_1.0.5         yaml_2.3.12             edgeR_4.8.2             snakecase_0.11.1       
#  [37] carData_3.0-6           qvalue_2.42.0           SparseArray_1.10.8      blob_1.3.0             
#  [41] crayon_1.5.3            ggtangle_0.1.1          lattice_0.22-7          cowplot_1.2.0          
#  [45] annotate_1.88.0         KEGGREST_1.50.0         pillar_1.11.1           knitr_1.51             
#  [49] fgsea_1.36.2            rjson_0.2.23            estimability_1.5.1      codetools_0.2-20       
#  [53] fastmatch_1.1-8         glue_1.8.0              ggiraph_0.9.6           ggfun_0.2.0            
#  [57] fontLiberation_0.1.0    vctrs_0.7.2             png_0.1-9               treeio_1.34.0          
#  [61] cellranger_1.1.0        gtable_0.3.6            cachem_1.1.0            xfun_0.55              
#  [65] openxlsx_4.2.8.1        S4Arrays_1.10.1         Seqinfo_1.0.0           coda_0.19-4.1          
#  [69] survival_3.8-3          iterators_1.0.14        statmod_1.5.1           ggtree_4.0.4           
#  [73] bit64_4.6.0-1           fontquiver_0.2.1        otel_0.2.0              colorspace_2.1-2       
#  [77] DBI_1.3.0               tidyselect_1.2.1        bit_4.6.0               compiler_4.5.1         
#  [81] httr2_1.2.2             fontBitstreamVera_0.1.1 DelayedArray_0.36.0     scales_1.4.0           
#  [85] rappdirs_0.3.4          digest_0.6.39           rmarkdown_2.31          XVector_0.50.0         
#  [89] htmltools_0.5.9         pkgconfig_2.0.3         fastmap_1.2.0           rlang_1.1.7            
#  [93] GlobalOptions_0.1.4     htmlwidgets_1.6.4       UCSC.utils_1.4.0        farver_2.1.2           
#  [97] jsonlite_2.0.0          GOSemSim_2.36.0         R.oo_1.27.1             magrittr_2.0.4         
# [101] Formula_1.2-5           GenomeInfoDbData_1.2.14 ggplotify_0.1.3         patchwork_1.3.2        
# [105] Rcpp_1.1.1              ape_5.8-1               ggnewscale_0.5.2        gdtools_0.5.0          
# [109] stringi_1.8.7           MASS_7.3-65             plyr_1.8.9              parallel_4.5.1         
# [113] listenv_0.10.1          Biostrings_2.78.0       splines_4.5.1           hms_1.1.4              
# [117] circlize_0.4.18         locfit_1.5-9.12         igraph_2.2.3            ggsignif_0.6.4         
# [121] reshape2_1.4.5          XML_3.99-0.23           evaluate_1.0.5          tzdb_0.5.0             
# [125] foreach_1.5.2           tweenr_2.0.3            polyclip_1.10-7         clue_0.3-68            
# [129] ggforce_0.5.0           xtable_1.8-8            tidytree_0.4.7          tidydr_0.0.6           
# [133] rstatix_0.7.3           aplot_0.2.9             memoise_2.0.1           AnnotationDbi_1.72.0   
# [137] cluster_2.1.8.2         timechange_0.4.0        globals_0.19.1         
