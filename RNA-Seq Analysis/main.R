library(recount)
library(DESeq2)
library(ggplot2)
library(dplyr)
library(tibble)
library(stringr)

library(EnhancedVolcano)
library(msigdbr)
library(EnsDb.Hsapiens.v86)
library(clusterProfiler)
library(enrichR)

# ---- Data ----

# Two groups: HeLa Ctrl and HeLa Salmonella infected MOI 100 (20 hours post-infection)
# Three biological replicates each

SELECTED_STUDY <- "SRP034009"
download_study(SELECTED_STUDY, outdir = "data")

RDATA_PATH <- file.path("data", "rse_gene.Rdata")

load(RDATA_PATH)

# ---- Question ----
# Salmonella infection least efficient in G1-arrested cells. Cyclin D1 promotes transition from G1 to S.
# miR-15 family of miRNAs decrease cyclin D1 activity.
# miR-15 miRNAs are downregulated on Salmonella infection.
# What genes are enriched by Salmonella infection and what miR-15 miRNA are downregulated?

# ---- Workflow ----
rse_gene$condition <- c(rep("HeLaCtrl", 3), 
                        rep("HeLaInfected", 3))
cnts <- assay(rse_gene)
exp_metadata <- as.data.frame(colData(rse_gene))

# split based on conditions and remove low-count genes
dds <- DESeqDataSet(rse_gene, design = ~condition)
dds <- dds[rowSums(counts(dds)) >= 25]

# https://www.biostars.org/p/445113/
# Use rlog to convert counts to regularized log scale for visualization purposes.
# PCA results based on features with highest variance.
# Linear scale would see mostly high mean features while log scale is inverse.
# Find balance based on library size.

# Without change to scale, small changes with low counts would dominate
# ex. A gene with 2 reads == 200% more reads than a gene with 1 reads
#     A gene with 101 reads == 1% more reads than a gene with 100 reads. 
rld <- rlog(dds, blind = TRUE)

# QC: PCA plot to see if clustering by some factor (could be group, sex, etc.)
# Condition explains variance in PC1
plot_pca_qc <- plotPCA(rld, intgroup = "condition")

ggsave(file.path("output", "plot_pca_qc.png"), plot = plot_pca_qc)

dds <- DESeq(dds)

# contrast is chr vector where 2nd and 3rd args are numerator and denominator.
# we want to see relative expression of Inf. vs. Ctrl
res <- results(dds, 
               contrast = c("condition", "HeLaInfected", "HeLaCtrl"), 
               parallel = TRUE)

# scatter plot of log 2 fold changes vs. mean of normalized counts
# How does fold change relate to how much gene is expressed.
# Variance at low cnts can be extreme
plot_ma_origfc <- DESeq2::plotMA(res)

# Shrink fold change - creates prior distribution of lfcs from genes with high counts. filters noise from low count genes.
# https://support.bioconductor.org/p/77461/#77466
# coef = 2 is the second item of resultsNames(dds)
resnorm <- lfcShrink(dds = dds, res = res, type = "normal", coef = 2)

plot_ma_shkfc <- DESeq2::plotMA(resnorm)

resdf <- as.data.frame(resnorm)
View(resdf)

# ---- Annotate Results ----
ens2sym <- AnnotationDbi::select(EnsDb.Hsapiens.v86, keys = keys(EnsDb.Hsapiens.v86),
                                 columns = c("SYMBOL"))

resdf <- resdf %>% 
  rownames_to_column(var = "ENSID") %>%
  mutate(GENEID = gsub(ENSID, pattern = "\\..+", replacement = "")) %>%
  dplyr::select(-ENSID) %>%
  inner_join(y = ens2sym, by = "GENEID")

# ---- Volcano Plot ----
plot_volcano <- EnhancedVolcano(resdf, lab = resdf$SYMBOL, 
                                pCutoff = 1e-100, FCcutoff = 3,
                                x = "log2FoldChange", y = "padj")

# ---- Heatmap of overexpressed and underexpressed ----

# ---- Over representation Analysis ----
# Not the best since only uses padj and logfoldchange to determine expression levels
# p cutoff is arbitrary and can be inflated based on total number of genes over/under expressed
over_expressed_genes <- resdf %>% 
  dplyr::filter(padj < 0.01 & log2FoldChange > 2) %>%
  pull(SYMBOL)

under_expressed_genes <- resdf %>% 
  dplyr::filter(padj < 0.01 & log2FoldChange < -2)%>%
  pull(SYMBOL)

# category C5 includes GO gene sets.
# gene set: gene
gene_sets <- msigdbr(species = "Homo sapiens", category = "C5") %>%
  dplyr::select(gs_name, gene_symbol)

egmt <- enricher(gene = over_expressed_genes, TERM2GENE = gene_sets)

edf <- as.data.frame(egmt)


# Gene Set Enrichment Analysis
resdf2 <- resdf %>%
  arrange(padj) %>%
  # handle inf when taking logs by choosing lowest possible value instead of 0.
  mutate(padj = case_when(padj == 0 ~ .Machine$double.xmin,
                          TRUE ~ padj)) %>%
  # calculate gsea metric
  mutate(gsea_metric = -log10(padj) * sign(log2FoldChange)) %>%
  dplyr::filter(!is.na(gsea_metric)) %>%
  # Order by GSEA
  arrange(desc(gsea_metric))
  
hist(resdf2$gsea_metric, breaks = 100)

ranks <- resdf2 %>%
  select(SYMBOL, gsea_metric) %>%
  distinct(SYMBOL, .keep_all = TRUE) %>%
  deframe()

gseas_res <- GSEA(geneList = ranks,
                  TERM2GENE = gene_sets)

gsea_res_df <- as.data.fram(gsea_res)


# Fig 3a
