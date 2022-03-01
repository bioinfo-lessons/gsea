rm(list = ls()) # R version 4.1.2
library(dplyr) # dplyr_1.0.7
library(stringr) # stringr_1.4.0

# --- Functions ---
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

output.cls <- function(metacolumn, filename) {
  out <- as.vector(metacolumn)
  lvls <- levels(factor(metacolumn, levels = unique(out)))
  cat(c(length(out), length(lvls), 1, "\n"), sep = " ", file = filename)
  cat(c("#", lvls, "\n"), sep = " ", file = filename, append = TRUE)
  cat(out, sep = " ", file = filename, append = TRUE)
}

# --- Data ---
normcounts <- read.table("input/normalizedcounts.tsv", header = TRUE, 
                         row.names = 1)

# --- Code ---
# Metadata column
meta <- colnames(normcounts)
meta <- str_remove(meta, "[0-9]*$") %>% # Remove trailing number
  str_replace(pattern = "m", replacement = "Male") %>%
  str_replace(pattern = "f", replacement = "Female") %>% # m = Male; f = Female
  factor(levels = c("Male", "Female")) # Make factor
meta

# Write GSEA input
output.gct(normcounts, filename = "gender.gct")
output.cls(meta, filename = "gender.cls")