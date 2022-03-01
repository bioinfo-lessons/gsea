rm(list = ls()) # R version 4.0.3
library(stringr) # stringr_1.4.0
library(DESeq2) # DESeq2_1.30.1
library(dplyr) # dplyr_1.0.7

# --- Functions ---
output.gmt <- function(l, filename) {
  if(file.exists(filename)) {
    warning(paste("Removing previous", filename, "file."))
    file.remove(filename)
  }
  invisible(lapply(1:length(l), function(i) {
    cat(c(names(l)[i], "na", l[[i]], "\n"), sep = "\t", file = filename,
        append = TRUE)
  }))
}

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
rnk <- data.frame(Feature = rownames(res.ape), LFC = res.ape$log2FoldChange)
head(rnk)
rnk$Feature <- str_remove(rnk$Feature, "\\..*$")
head(rnk)

# Save .rnk (without header and tab separated)
write.table(rnk, file = "airway.rnk", sep = "\t", quote = FALSE, 
            col.names = FALSE, row.names = FALSE)

# How to create a gmt with my DEA results? Pick the top/bottom 100 (or 150, 
# 200,...) ordered by LFC (independently of the p-value)
up100 <- rnk %>% arrange(desc(LFC))
head(up100)
up100 <- up100$Feature[1:100]

down100 <- rnk %>% arrange(LFC)
head(down100)
down100 <- down100$Feature[1:100]

geneset <- list(up_airway = up100, down_airway = down100)
output.gmt(geneset, filename = "airway.gmt")