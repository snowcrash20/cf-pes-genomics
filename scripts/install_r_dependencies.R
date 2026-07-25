#!/usr/bin/env Rscript

cran_packages <- c(
  "dplyr",
  "ggplot2",
  "purrr",
  "readr",
  "rlang",
  "stringr",
  "tibble",
  "tidyr"
)

bioconductor_packages <- c(
  "Biostrings",
  "GenomicRanges",
  "IRanges",
  "rtracklayer",
  "S4Vectors"
)

missing_cran <- cran_packages[!vapply(
  cran_packages,
  requireNamespace,
  quietly = TRUE,
  FUN.VALUE = logical(1)
)]

if (length(missing_cran) > 0L) {
  install.packages(missing_cran, repos = "https://cloud.r-project.org")
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

missing_bioconductor <- bioconductor_packages[!vapply(
  bioconductor_packages,
  requireNamespace,
  quietly = TRUE,
  FUN.VALUE = logical(1)
)]

if (length(missing_bioconductor) > 0L) {
  BiocManager::install(missing_bioconductor, ask = FALSE, update = FALSE)
}

message("R dependencies are installed.")
