#!/usr/bin/env Rscript

annotated <- readr::read_csv(
  "results/example/snp/tables/annotated_snps.csv",
  show_col_types = FALSE
)
filter_summary <- readr::read_csv(
  "results/example/snp/tables/filter_summary.csv",
  show_col_types = FALSE
)
snp_status <- readr::read_csv(
  "results/example/snp/tables/snp_status.csv",
  show_col_types = FALSE
)
contigs <- readr::read_csv(
  "results/example/contigs/tables/contig_taxonomy.csv",
  show_col_types = FALSE
)

stopifnot(
  nrow(annotated) == 10L,
  sum(annotated$feature_label == "Intergenic") == 2L,
  identical(filter_summary$retained_snps, c(4, 5)),
  identical(snp_status$Fixed, c(1, 1)),
  identical(snp_status$Marginal, c(3, 4)),
  nrow(contigs) == 4L,
  sum(contigs$length_bp) == 370,
  sum(contigs$taxon_name == "Unclassified") == 1L
)

message("Synthetic output values are correct.")
