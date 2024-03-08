rm(list = ls()) # R version 4.3.1 (2023-06-16)
library(tidyverse) # tidyverse_2.0.0

# --- Functions ---
# Creates a GCT file to load a normalized expression matrix to GSEA
output.gct <- function(ncounts, filename) {
  NAME <- Description <- rownames(ncounts)
  out <- cbind(NAME, Description, ncounts)
  cat(c("#1.2", "\n"), sep = "\t", file = filename)
  cat(c(nrow(ncounts), ncol(ncounts), "\n"), sep = "\t", file = filename,
      append = TRUE)
  suppressWarnings(
    write.table(out, file = filename, row.names = FALSE, quote = FALSE,
                sep = "\t", na = "", append = TRUE))
}

# Creates a CLS file to load the phenotype labels to GSEA
output.cls <- function(metacolumn, filename) {
  out <- as.vector(metacolumn)
  lvls <- levels(factor(metacolumn, levels = unique(out)))
  cat(c(length(out), length(lvls), 1, "\n"), sep = " ", file = filename)
  cat(c("#", lvls, "\n"), sep = " ", file = filename, append = TRUE)
  cat(out, sep = " ", file = filename, append = TRUE)
}

# --- Data ---
# Normalized expression counts
normcounts <- read.table("input/normalizedcounts.tsv", header = TRUE,
                         row.names = 1)

# --- Code ---
# Create a metadata factor
meta <- colnames(normcounts) # Sample names
meta <- str_remove(meta, "[0-9]*$") %>% # Remove trailing number
  str_replace(pattern = "m", replacement = "Male") %>%
  str_replace(pattern = "f", replacement = "Female") %>% # m = Male; f = Female
  factor(levels = c("Male", "Female")) # Create a factor
meta

# Write GSEA input
output.gct(normcounts, filename = "gender.gct")
output.cls(meta, filename = "gender.cls")
