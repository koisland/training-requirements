# ----- Phantasus RNASeq Visulation in Browser -----
# library(phantasus)

# servePhantasus(host = "127.0.0.1")

# ----- Examples: Airway and DESeq -----
library(airway)
library(DESeq2)
library(ggplot2)
library(dplyr)

data("airway")
se <- airway

thr_params <- list(pThr = 0.01, logFCThr = 1, baseMeanThr = 20, cpmThr = 1)
thr_params_df <- as.data.frame(thr_params)

head(assay(se))
rowData(se)
colData(se)

# pairwise comparison of cell types with treatment
dds <- DESeqDataSet(se, design = ~ cell + dex)

# run diff. gene expression analysis
dds <- DESeq(dds)

res <- results(dds)

# select significant genes.
sigRes <- as.data.frame(res) %>%
  filter(padj <= thr_params$pThr) %>%
  filter(abs(log2FoldChange) >= thr_params$logFCThr) %>%
  filter(baseMean >= thr_params$baseMeanThr)

write.table(thr_params_df, file = "output/airway_sigGenes_thrs_1.csv",  sep = ",")
write.table(sigRes, file = "output/airway_sigGenes_1.csv", sep = ",", col.names = TRUE, row.names = TRUE)
