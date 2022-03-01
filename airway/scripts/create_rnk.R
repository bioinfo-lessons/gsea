rm(list = ls()) # R version 4.1.2
library(stringr) # stringr_1.4.0
library(DESeq2) # DESeq2_1.34.0
library(dplyr)

# --- Functions ---
output.gmt <- function(l, filename) {
  if(file.exists(filename)) {
    warning(paste("Removing previous", filename, "file."))
    file.remove(filename)
  }
  invisible(lapply(seq_along(l), function(i) {
    cat(c(names(l)[i], "na", l[[i]], "\n"), sep = "\t", file = filename,
        append = TRUE)
  }))
}

# --- Data ---
dds <- readRDS("input/dds.rds")

# --- Code ---
# DEGs
res <- results(dds, alpha = 0.05, contrast = c("dex", "trt", "untrt"))
summary(res)
res

# Shrunken LFC
resultsNames(dds)
res.ape <- lfcShrink(dds, coef = "dex_trt_vs_untrt", type = "apeglm",
                     res = res)
head(res.ape)

# Create .rnk
rnk <- data.frame(Feature = rownames(res.ape), LFC = res.ape$log2FoldChange)
head(rnk)
rnk$Feature <- str_remove(rnk$Feature, "\\..*$")

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