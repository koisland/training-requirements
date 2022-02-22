library(recount)
library(DESeq2)
library(ggplot2)
library(dplyr)

library(EnhancedVolcano)
library(msigdbr)

# ---- Get Data ----
selected_study <- recount::abstract_search("Pneumolysin")$project

download_study(selected_study)

rdata_path <- file.path(selected_study, "rse_gene.Rdata")

load(rdata_path)

# ---- Original workflow (Human transcriptomic analysis) ----
# Align reads in paired-end mode to human genome using STAR aligner
# Mapped data (BAM?) converted to gene counts using HTSeqcount and UCSC annotation
# DESeq for DEA at every time point comparing infected to untreated samples

# ---- Question ----
# How do human airway epithelial cells respond to statin/cytolysin challenge?
# Four different classes from three different patients: vehicle control, only statin, only toxin, and both statin and toxin

# ---- Workflow ----
sample_types <- c(rep("EtOH_PBS_Control", 3), 
                  rep("SimvastatinOnly", 2), 
                  rep("BOTHSimvastatinPLY", 3),
                  rep("PLY_Only", 3),
                  "SimvastatinOnly")

rse_gene$condition <- sample_types


cnts <- assay(rse_gene)
exp_metadata <- as.data.frame(colData(rse_gene))

# split based on condition
dds <- DESeqDataSet(rse_gene, design = ~condition)

# https://www.biostars.org/p/445113/
# Use rlog to convert counts to log scale.
# PCA results based on features with highest variance.
# Linear scale would see mostly high mean features while log scale is inverse.

# Without change to scale, small changes with low counts would dominate
# ex. A gene with 2 reads == 200% more reads than a gene with 1 reads
#     A gene with 101 reads == 1% more reads than a gene with 100 reads. 
rld <- rlog(dds)

plotPCA(rld)

dds <- DESeq(dds)

# contrast is chr vector where 2nd and 3rd args are numerator and denominator.
# we want to see relative expression of ctr vs. inf
res <- results(dds, contrast = c("condition", ))

# scatter plot of log 2 fold changes vs. mean of normalized counts
# 
plotMA(res)

# Change times where fold change
resnorm <- lfcShrink(dds = dds, res = res, type = "normal", coef = 2)

plotMA(resnorm)

resdf <- as.data.frame(resnorm)
View(resdf)
# ---- Annotate Report ----