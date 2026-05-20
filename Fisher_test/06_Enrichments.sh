# Enrichment.r
# -g = Two columns table of input genes with specific association from your study
# -l = list of two columns tables with gene - disease association. E.g. Gene1 SYN
# -p = make a bubble chart with OR and -log10(FDR)
# -b = background (protein coding = 19776, brain expressed = 15585, WGCNA list = 6029)
# -o = output label for statistics and viz
# -W/-H = width/height of the plot. 

mkdir enrichments/
cp utils/geneset/*.RData enrichments/
cp utils/Enrichment.r enrichments/
cp deg.Combo.DMSO.txt enrichments/
cp deg.BMS.DMSO.txt enrichments/
cp deg.Ribo.DMSO.txt enrichments/

cd enrichments/

mkdir STATS/

## Running Script for Combo, Ribo and BMS vs DMSO

Rscript Enrichment.r -g deg.Combo.DMSO.txt -l GeneSets_TME_Palbociclib_DMSO.RData -p -b 15585 -o STATS/Pritika_project_TME_markers_Combo -W 5 -H 6
Rscript Enrichment.r -g deg.BMS.DMSO.txt -l GeneSets_TME_Palbociclib_DMSO.RData -p -b 15585 -o STATS/Pritika_project_TME_markers_BMS -W 5 -H 6
Rscript Enrichment.r -g deg.Ribo.DMSO.txt -l GeneSets_TME_Palbociclib_DMSO.RData -p -b 15585 -o STATS/Pritika_project_TME_markers_Ribo -W 5 -H 6


rm *.RData
