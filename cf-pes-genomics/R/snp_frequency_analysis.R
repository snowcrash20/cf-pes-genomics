#!/usr/bin/env Rscript

usage <- paste(
  "SNP frequency analysis for breseq-derived SNP tables",
  "",
  "Usage:",
  "  Rscript R/snp_frequency_analysis.R \\",
  "    --snp-dir PATH \\",
  "    --gff PATH \\",
  "    --metadata PATH \\",
  "    --out-dir PATH \\",
  "    [--aliases PATH] \\",
  "    [--qual-threshold 10] \\",
  "    [--fixed-threshold 0.95] \\",
  "    [--top-n 5]",
  "",
  "Input SNP files must be named <sample>_snps_breseq.tsv and contain:",
  "CHROM, POS, REF, ALT, QUAL, AF, AD, DP.",
  sep = "\n"
)

raw_args <- commandArgs(trailingOnly = TRUE)

if (any(raw_args %in% c("-h", "--help"))) {
  cat(usage, "\n")
  quit(status = 0L)
}

file_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(file_argument) == 0L) {
  normalizePath("R", mustWork = TRUE)
} else {
  dirname(normalizePath(sub("^--file=", "", file_argument[[1]]), mustWork = TRUE))
}
source(file.path(script_dir, "common.R"))

options <- parse_cli_args(
  raw_args,
  defaults = list(
    aliases = NULL,
    qual_threshold = "10",
    fixed_threshold = "0.95",
    top_n = "5"
  )
)

require_cli_values(options, c("snp_dir", "gff", "metadata", "out_dir"))
require_packages(c(
  "dplyr",
  "GenomicRanges",
  "ggplot2",
  "IRanges",
  "purrr",
  "readr",
  "rlang",
  "rtracklayer",
  "S4Vectors",
  "stringr",
  "tibble",
  "tidyr"
))

qual_threshold <- parse_number(
  options$qual_threshold,
  "qual_threshold",
  minimum = 0
)
fixed_threshold <- parse_number(
  options$fixed_threshold,
  "fixed_threshold",
  minimum = 0,
  maximum = 1
)
top_n <- parse_integer(options$top_n, "top_n")

require_directory(options$snp_dir, "SNP directory")
require_file(options$gff, "GFF3 annotation")

metadata <- read_sample_metadata(options$metadata)

dir.create(options$out_dir, recursive = TRUE, showWarnings = FALSE)
table_dir <- file.path(options$out_dir, "tables")
figure_dir <- file.path(options$out_dir, "figures")
comparison_dir <- file.path(figure_dir, "af_comparisons")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(comparison_dir, recursive = TRUE, showWarnings = FALSE)

read_snp_table <- function(path) {
  raw <- readr::read_tsv(
    path,
    col_names = FALSE,
    col_types = readr::cols(.default = readr::col_character()),
    show_col_types = FALSE,
    progress = FALSE,
    comment = ""
  )

  if (ncol(raw) != 8L) {
    stop(
      basename(path),
      " has ",
      ncol(raw),
      " columns; expected exactly 8.",
      call. = FALSE
    )
  }

  if (nrow(raw) > 0L && toupper(sub("^#", "", raw[[1]][[1]])) == "CHROM") {
    raw <- dplyr::slice(raw, -1L)
  }

  names(raw) <- c("CHROM", "POS", "REF", "ALT", "QUAL", "AF", "AD", "DP")
  sample_id <- stringr::str_remove(
    basename(path),
    "_snps_breseq\\.tsv$"
  )

  parsed <- dplyr::mutate(
    raw,
    CHROM = as.character(CHROM),
    POS = suppressWarnings(as.integer(POS)),
    REF = as.character(REF),
    ALT = as.character(ALT),
    QUAL = suppressWarnings(as.numeric(QUAL)),
    AF = suppressWarnings(as.numeric(AF)),
    AD = as.character(AD),
    DP = suppressWarnings(as.numeric(DP)),
    sample = sample_id
  )

  invalid <- is.na(parsed$POS) | is.na(parsed$QUAL) | is.na(parsed$AF)
  if (any(invalid)) {
    stop(
      basename(path),
      " contains non-numeric POS, QUAL, or AF values.",
      call. = FALSE
    )
  }

  if (any(parsed$AF < 0 | parsed$AF > 1)) {
    stop(
      basename(path),
      " contains AF values outside the interval [0, 1].",
      call. = FALSE
    )
  }

  parsed
}

snp_files <- list.files(
  options$snp_dir,
  pattern = "_snps_breseq\\.tsv$",
  full.names = TRUE
)

if (length(snp_files) == 0L) {
  stop(
    "No files matching '*_snps_breseq.tsv' were found in ",
    options$snp_dir,
    ".",
    call. = FALSE
  )
}

all_snps <- purrr::map_dfr(snp_files, read_snp_table) |>
  dplyr::mutate(snp_row_id = dplyr::row_number())

missing_metadata <- setdiff(unique(all_snps$sample), metadata$sample)
if (length(missing_metadata) > 0L) {
  stop(
    "These SNP samples are absent from the metadata file: ",
    paste(missing_metadata, collapse = ", "),
    call. = FALSE
  )
}

missing_snp_files <- setdiff(metadata$sample, unique(all_snps$sample))
if (length(missing_snp_files) > 0L) {
  stop(
    "These metadata samples have no matching SNP file: ",
    paste(missing_snp_files, collapse = ", "),
    call. = FALSE
  )
}

gff <- rtracklayer::import(options$gff)
feature_type <- as.character(S4Vectors::mcols(gff)$type)
keep <- feature_type %in% c("CDS", "gene")

if (!any(keep)) {
  stop("The GFF3 contains no CDS or gene features.", call. = FALSE)
}

features <- gff[keep]
feature_type <- feature_type[keep]

metadata_as_character <- function(ranges, candidates) {
  result <- rep(NA_character_, length(ranges))
  columns <- colnames(S4Vectors::mcols(ranges))

  for (candidate in candidates[candidates %in% columns]) {
    values <- S4Vectors::mcols(ranges)[[candidate]]
    values <- vapply(
      as.list(values),
      function(value) paste(as.character(value), collapse = ";"),
      character(1)
    )
    values[!nzchar(values)] <- NA_character_
    replace <- is.na(result) & !is.na(values)
    result[replace] <- values[replace]
  }

  result
}

feature_annotations <- tibble::tibble(
  feature_index = seq_along(features),
  feature_type = feature_type,
  gene = metadata_as_character(
    features,
    c("gene", "Name", "locus_tag", "ID")
  ),
  product = metadata_as_character(features, c("product", "description"))
)

snp_ranges <- GenomicRanges::GRanges(
  seqnames = all_snps$CHROM,
  ranges = IRanges::IRanges(
    start = all_snps$POS,
    width = pmax(nchar(all_snps$REF), 1L)
  )
)

hits <- GenomicRanges::findOverlaps(
  snp_ranges,
  features,
  ignore.strand = TRUE
)

if (length(hits) > 0L) {
  chosen_annotations <- tibble::tibble(
    snp_row_id = S4Vectors::queryHits(hits),
    feature_index = S4Vectors::subjectHits(hits)
  ) |>
    dplyr::left_join(feature_annotations, by = "feature_index") |>
    dplyr::mutate(feature_priority = dplyr::if_else(feature_type == "CDS", 1L, 2L)) |>
    dplyr::arrange(snp_row_id, feature_priority, feature_index) |>
    dplyr::group_by(snp_row_id) |>
    dplyr::slice_head(n = 1L) |>
    dplyr::ungroup() |>
    dplyr::select(snp_row_id, feature_type, gene, product)
} else {
  chosen_annotations <- tibble::tibble(
    snp_row_id = integer(),
    feature_type = character(),
    gene = character(),
    product = character()
  )
}

annotated_snps <- all_snps |>
  dplyr::left_join(chosen_annotations, by = "snp_row_id")

if (!is.null(options$aliases)) {
  require_file(options$aliases, "Product alias file")
  aliases <- readr::read_tsv(
    options$aliases,
    show_col_types = FALSE,
    progress = FALSE,
    col_types = readr::cols(.default = readr::col_character())
  )

  if (!all(c("product", "label") %in% names(aliases))) {
    stop("Alias file must contain 'product' and 'label' columns.", call. = FALSE)
  }

  aliases <- dplyr::distinct(aliases, product, .keep_all = TRUE)
  annotated_snps <- dplyr::left_join(
    annotated_snps,
    dplyr::rename(aliases, product_alias = label),
    by = "product"
  )
} else {
  annotated_snps$product_alias <- NA_character_
}

annotated_snps <- annotated_snps |>
  dplyr::mutate(
    feature_label = dplyr::coalesce(
      product_alias,
      gene,
      product,
      "Intergenic"
    ),
    snp_id = paste(CHROM, POS, REF, ALT, sep = "_")
  ) |>
  dplyr::left_join(metadata, by = "sample") |>
  dplyr::arrange(order, CHROM, POS)

filtered_snps <- dplyr::filter(annotated_snps, QUAL >= qual_threshold)

filter_summary <- annotated_snps |>
  dplyr::group_by(sample, label, patient, timepoint, order) |>
  dplyr::summarise(
    total_snps = dplyr::n(),
    retained_snps = sum(QUAL >= qual_threshold, na.rm = TRUE),
    filtered_out = total_snps - retained_snps,
    percent_removed = dplyr::if_else(
      total_snps > 0L,
      100 * filtered_out / total_snps,
      0
    ),
    .groups = "drop"
  ) |>
  dplyr::arrange(order)

snp_status_counts <- filtered_snps |>
  dplyr::mutate(
    status = dplyr::if_else(AF >= fixed_threshold, "Fixed", "Marginal")
  ) |>
  dplyr::count(
    sample,
    label,
    patient,
    timepoint,
    order,
    status,
    name = "count"
  )

snp_status <- tidyr::crossing(
  dplyr::select(metadata, sample, label, patient, timepoint, order),
  status = c("Fixed", "Marginal")
) |>
  dplyr::left_join(
    snp_status_counts,
    by = c("sample", "label", "patient", "timepoint", "order", "status")
  ) |>
  dplyr::mutate(count = tidyr::replace_na(count, 0L)) |>
  tidyr::pivot_wider(
    names_from = status,
    values_from = count,
    values_fill = 0L
  ) |>
  dplyr::mutate(Total = Fixed + Marginal) |>
  dplyr::arrange(order)

top_features <- filtered_snps |>
  dplyr::filter(feature_label != "Intergenic") |>
  dplyr::count(
    sample,
    label,
    order,
    feature_label,
    name = "snp_count"
  ) |>
  dplyr::group_by(sample) |>
  dplyr::slice_max(snp_count, n = top_n, with_ties = FALSE) |>
  dplyr::ungroup() |>
  dplyr::arrange(order, dplyr::desc(snp_count))

readr::write_csv(
  dplyr::select(annotated_snps, -snp_row_id),
  file.path(table_dir, "annotated_snps.csv")
)
readr::write_csv(filter_summary, file.path(table_dir, "filter_summary.csv"))
readr::write_csv(snp_status, file.path(table_dir, "snp_status.csv"))
readr::write_csv(top_features, file.path(table_dir, "top_features.csv"))

if (nrow(top_features) > 0L) {
  ordered_labels <- metadata |>
    dplyr::arrange(order) |>
    dplyr::pull(label)

  top_features <- dplyr::mutate(
    top_features,
    label = factor(label, levels = ordered_labels)
  )

  top_plot <- ggplot2::ggplot(
    top_features,
    ggplot2::aes(x = label, y = snp_count, fill = feature_label)
  ) +
    ggplot2::geom_col() +
    ggplot2::labs(
      title = paste0("Top ", top_n, " annotated features per sample"),
      subtitle = paste0("SNPs retained at QUAL >= ", qual_threshold),
      x = "Sample",
      y = "SNP count",
      fill = "Feature"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid.major.x = ggplot2::element_blank()
    )

  ggplot2::ggsave(
    file.path(figure_dir, "top_features_per_sample.png"),
    top_plot,
    width = 11,
    height = 7,
    dpi = 300
  )
}

make_pair_table <- function(patient_metadata) {
  patient_metadata <- dplyr::arrange(patient_metadata, order)

  if (nrow(patient_metadata) < 2L) {
    return(tibble::tibble())
  }

  indices <- utils::combn(seq_len(nrow(patient_metadata)), 2L)
  tibble::tibble(
    patient = patient_metadata$patient[indices[1L, ]],
    sample_x = patient_metadata$sample[indices[1L, ]],
    sample_y = patient_metadata$sample[indices[2L, ]],
    label_x = patient_metadata$timepoint[indices[1L, ]],
    label_y = patient_metadata$timepoint[indices[2L, ]]
  )
}

comparison_pairs <- metadata |>
  split(metadata$patient) |>
  purrr::map_dfr(make_pair_table)

plot_af_comparison <- function(patient, sample_x, sample_y, label_x, label_y) {
  pair_snps_long <- filtered_snps |>
    dplyr::filter(sample %in% c(sample_x, sample_y)) |>
    dplyr::group_by(snp_id, sample, feature_label) |>
    dplyr::summarise(AF = mean(AF), .groups = "drop")

  if (nrow(pair_snps_long) == 0L) {
    warning(
      "Skipping ",
      patient,
      " ",
      label_x,
      " versus ",
      label_y,
      ": neither sample retained a SNP."
    )
    return(invisible(NULL))
  }

  pair_snps <- pair_snps_long |>
    tidyr::pivot_wider(
      names_from = sample,
      values_from = AF,
      values_fill = 0
    )

  for (sample_name in c(sample_x, sample_y)) {
    if (!sample_name %in% names(pair_snps)) {
      pair_snps[[sample_name]] <- 0
    }
  }

  pair_snps <- dplyr::rename(
    pair_snps,
    AF_x = !!rlang::sym(sample_x),
    AF_y = !!rlang::sym(sample_y)
  )

  top_labels <- filtered_snps |>
    dplyr::filter(
      sample %in% c(sample_x, sample_y),
      feature_label != "Intergenic"
    ) |>
    dplyr::count(feature_label, name = "count") |>
    dplyr::slice_max(count, n = top_n, with_ties = FALSE) |>
    dplyr::pull(feature_label)

  pair_snps <- dplyr::mutate(
    pair_snps,
    colour_group = dplyr::if_else(
      feature_label %in% top_labels,
      feature_label,
      "Other"
    ),
    point_alpha = dplyr::if_else(colour_group == "Other", 0.35, 0.9)
  )

  comparison_plot <- ggplot2::ggplot(
    pair_snps,
    ggplot2::aes(x = AF_x, y = AF_y)
  ) +
    ggplot2::geom_point(
      ggplot2::aes(color = colour_group, alpha = point_alpha),
      size = 1.4
    ) +
    ggplot2::scale_alpha_identity() +
    ggplot2::geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed",
      color = "grey45"
    ) +
    ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    ggplot2::labs(
      title = patient,
      subtitle = paste(label_x, "versus", label_y),
      x = paste0(label_x, " allele frequency"),
      y = paste0(label_y, " allele frequency"),
      color = "Feature"
    ) +
    ggplot2::theme_minimal(base_size = 12)

  filename <- paste0(
    slugify(patient),
    "_",
    slugify(label_x),
    "_vs_",
    slugify(label_y),
    ".png"
  )

  ggplot2::ggsave(
    file.path(comparison_dir, filename),
    comparison_plot,
    width = 7,
    height = 6,
    dpi = 300
  )
}

if (nrow(comparison_pairs) > 0L) {
  purrr::pwalk(comparison_pairs, plot_af_comparison)
}

write_session_info(file.path(options$out_dir, "sessionInfo.txt"))

message(
  "Processed ",
  length(snp_files),
  " SNP files. Outputs: ",
  normalizePath(options$out_dir, mustWork = TRUE)
)
