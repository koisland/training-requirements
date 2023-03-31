library(dplyr)
library(ggplot2)

p_thr <- 0.01
log_fold_change_thr <- 2

rse <- readRDS("EwS.rds")

# Sample metadata.
# We have column condition to split our dataset by.
cd <- as.data.frame(
  SummarizedExperiment::colData(rse)
)
testit::assert("condition" %in% colnames(cd))

# Design is what column to use to split our dataset
deseq_dataset <- DESeq2::DESeqDataSet(rse, design = ~condition)

# (1)
# Log transform count data.
# Plot PCA of counts.
# First PC (PC1) cleanly splits samples by type.
# * Control vs. EF1 (EWSR1-FLI1) knockdown.
pca <- DESeq2::rlog(deseq_dataset) %>% BiocGenerics::plotPCA()

# Based on PCA, we have good data
# * Treatment not artifacts (batch effects) responsible for technical variation.
# Perform DE analysis.
dds <- DESeq2::DESeq(deseq_dataset)

# contrast is chr vector where 2nd and 3rd args are numerator and denominator.
# We want to see relative expression of Knockdown. vs. Ctrl
res <- DESeq2::results(dds, contrast = c("condition", "shEF1", "shCTR"))

# scatter plot of log 2 fold changes vs. mean of normalized counts
# How does fold change relate to how much gene is expressed.
# Variance at low cnts can be extreme. See upper left corner.
plot_ma_origfc <- DESeq2::plotMA(res)

# Shrink fold change -
# * creates prior distribution of lfcs from genes with high counts.
# * filters noise from low count genes.
# https://support.bioconductor.org/p/77461/#77466
# coef = 2 is the second item of resultsNames(dds)
res_shrunkfc <- DESeq2::lfcShrink(
  dds = dds,
  res = res,
  type = "normal",
  coef = 2
)

# (2)
# Then replot with log fc shrunk results.
plot_ma_shrunkfc <- DESeq2::plotMA(res_shrunkfc)

# Annotate deseq gene IDs with symbols.
annotate_deseq_res <- function(resdf) {
  ens2sym <- AnnotationDbi::select(
    EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86,
    keys = AnnotationDbi::keys(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86),
    columns = c("SYMBOL")
  )

  return(
    # nolint start
    resdf %>%
      tibble::rownames_to_column(var = "ENSID") %>%
      # Remove any trailing digits from IDs. ex. ########.##
      dplyr::mutate(GENEID = gsub(ENSID, pattern = "\\..+", replacement = "")) %>%
      dplyr::select(-ENSID) %>%
      dplyr::inner_join(y = ens2sym, by = "GENEID")
    # nolint end
  )
}

resdf <- as.data.frame(res_shrunkfc) %>% annotate_deseq_res()

# (4)
# Visual identification of genes with large fold changes that are significant.
plot_volcano <- EnhancedVolcano::EnhancedVolcano(
  resdf,
  lab = resdf$SYMBOL,
  pCutoff = 1e-100,
  FCcutoff = 3,
  x = "log2FoldChange", y = "padj"
)

# (3)
resdf_de_genes <- resdf %>%
  # Must be significant.
  dplyr::filter(padj < p_thr) %>%
  # Add column to describe expression level.
  dplyr::mutate(
    type = dplyr::case_when(
      log2FoldChange > log_fold_change_thr ~ "overexpressed",
      log2FoldChange < -log_fold_change_thr ~ "underexpressed",
      .default = "other"
    )
  ) %>%
  # Sort by log2 fold change.
  dplyr::arrange(desc(log2FoldChange)) %>%
  # And add a rank.
  dplyr::mutate(rank = row_number())

num_de_genes <- nrow(resdf_de_genes)
resdf_top_de_genes <- resdf_de_genes %>%
  dplyr::filter(
    between(rank, 1, 10) | between(rank, num_de_genes - 9, num_de_genes)
  )

# (5)
plot_top_de_genes <- ggplot2::ggplot(
  resdf_top_de_genes,
  ggplot2::aes(
    x = log2FoldChange,
    y = SYMBOL,
    fill = type,
    label = sprintf("%0.2f", round(log2FoldChange, digits = 2))
  )
) +
  ggplot2::ggtitle("Top 10 Differentially Expressed Genes") +
  ggplot2::geom_bar(stat = "identity") +
  ggplot2::geom_text(vjust = "center") +
  ggplot2::theme_classic()

# Get over and under expressed genes.
overexp_genes <- resdf_de_genes %>%
  dplyr::filter(type == "overexpressed") %>%
  dplyr::pull(SYMBOL)

underexp_genes <- resdf_de_genes %>%
  dplyr::filter(type == "underexpressed") %>%
  dplyr::pull(SYMBOL)

# Get gene sets from H. sapiens
# Select the gene set name and the gene symbol.
gene_sets <- msigdbr::msigdbr(species = "Homo sapiens", category = "C5") %>%
  dplyr::select(gs_name, gene_symbol)

egmt <- clusterProfiler::enricher(
  gene = overexp_genes,
  TERM2GENE = gene_sets
)

edf <- as.data.frame(egmt)

clusterProfiler::dotplot(egmt)
barplot(egmt)
