rm(list = ls()) # R version 4.3.1 (2023-06-16)
library(tidyverse) # tidyverse_2.0.0
library(DESeq2) # DESeq2_1.42.0

# --- Data ---
dds <- readRDS("input/dds.rds")

# --- Code ---
# Review the metadata
colData(dds)
dds$cell
dds$dex

# Review the design
design(dds)

# DEGs
res <- results(dds, alpha = 0.05, contrast = c("dex", "trt", "untrt"))
summary(res)
res

#If you decide use LFC
# Shrunken LFC: Necessary for ranking!!! Also useful to compare LFC across 
# experiments. It is recommended to use these shrunken results when doing DEA.
# Please note that the p-values might vary a little.

### Low count genes have a lot of variance and thus their estimated LFCs can be 
### very extreme and imprecise.
### lfcShrink allows to shrunk the LFC of genes with low counts while conserving
### the LFC of genes with a high number of counts.
### apeglm is the default method, but there are others available:
### https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#moreshrink
resultsNames(dds)
res.ape <- lfcShrink(dds, coef = "dex_trt_vs_untrt", type = "apeglm",
                     res = res)
summary(res.ape) # Same number of up/down genes padj < 0.05

# Visual explanation of shrunken LFCs
par(mfrow = c(1, 2))
plotMA(res, ylim = c(-3, 3))
plotMA(res.ape, ylim = c(-3, 3))

# Create .rnk
rnk <- data.frame(Feature = rownames(res), Stat = res$stat)
head(rnk)
rnk$Feature <- str_remove(rnk$Feature, "\\..*$")
head(rnk)

# Save .rnk (without header and tab separated)
write.table(rnk, file = "airway.rnk", sep = "\t", quote = FALSE, 
            col.names = FALSE, row.names = FALSE)
