#!/usr/bin/env Rscript

usage <- paste(
  "Join MEGAHIT contigs to Kraken2 classifications and summarise their origin",
  "",
  "Usage:",
  "  Rscript R/unmapped_contig_taxonomy.R \\",
  "    --megahit-dir PATH \\",
  "    --kraken-dir PATH \\",
  "    --metadata PATH \\",
  "    --out-dir PATH \\",
  "    [--top-n 7]",
  "",
  "Expected files:",
  "  <megahit-dir>/<sample>/final.contigs.fa",
  "  <kraken-dir>/<sample>.out",
  "  <kraken-dir>/<sample>.report",
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
  defaults = list(top_n = "7")
)

require_cli_values(
  options,
  c("megahit_dir", "kraken_dir", "metadata", "out_dir")
)
require_packages(c(
  "Biostrings",
  "dplyr",
  "ggplot2",
  "purrr",
  "readr",
  "stringr",
  "tibble",
  "tidyr"
))

top_n <- parse_integer(options$top_n, "top_n")
require_directory(options$megahit_dir, "MEGAHIT directory")
require_directory(options$kraken_dir, "Kraken2 directory")

metadata <- read_sample_metadata(options$metadata) |>
  dplyr::arrange(order)

dir.create(options$out_dir, recursive = TRUE, showWarnings = FALSE)
table_dir <- file.path(options$out_dir, "tables")
figure_dir <- file.path(options$out_dir, "figures")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

read_kraken_report <- function(path) {
  raw <- readr::read_tsv(
    path,
    col_names = FALSE,
    col_types = readr::cols(.default = readr::col_character()),
    show_col_types = FALSE,
    progress = FALSE,
    trim_ws = FALSE
  )

  if (ncol(raw) == 6L) {
    names(raw) <- c(
      "pct",
      "clade_reads",
      "taxon_reads",
      "rank",
      "taxid",
      "taxon_name"
    )
  } else if (ncol(raw) == 8L) {
    names(raw) <- c(
      "pct",
      "clade_reads",
      "taxon_reads",
      "minimizers",
      "distinct_minimizers",
      "rank",
      "taxid",
      "taxon_name"
    )
  } else {
    stop(
      basename(path),
      " has ",
      ncol(raw),
      " columns; expected a standard 6-column or minimizer 8-column report.",
      call. = FALSE
    )
  }

  raw |>
    dplyr::transmute(
      taxid = as.character(taxid),
      taxon_name_report = stringr::str_trim(taxon_name),
      rank = as.character(rank)
    ) |>
    dplyr::distinct(taxid, .keep_all = TRUE)
}

read_kraken_output <- function(path) {
  raw <- readr::read_tsv(
    path,
    col_names = FALSE,
    col_types = readr::cols(.default = readr::col_character()),
    show_col_types = FALSE,
    progress = FALSE,
    quote = "",
    comment = ""
  )

  if (ncol(raw) != 5L) {
    stop(
      basename(path),
      " has ",
      ncol(raw),
      " columns; expected standard 5-column Kraken2 output.",
      call. = FALSE
    )
  }

  names(raw) <- c(
    "classification",
    "contig_id",
    "taxon_field",
    "kraken_length",
    "hit_list"
  )

  named_taxon <- stringr::str_detect(
    raw$taxon_field,
    "\\(taxid\\s+[0-9]+\\)$"
  )

  raw |>
    dplyr::transmute(
      classification = as.character(classification),
      contig_id = as.character(contig_id),
      taxid = dplyr::if_else(
        named_taxon,
        stringr::str_match(taxon_field, "taxid\\s+([0-9]+)")[, 2],
        as.character(taxon_field)
      ),
      taxon_name_inline = dplyr::if_else(
        named_taxon,
        stringr::str_trim(
          stringr::str_remove(taxon_field, "\\s*\\(taxid.*$")
        ),
        NA_character_
      ),
      kraken_length = suppressWarnings(as.numeric(kraken_length))
    ) |>
    dplyr::distinct(contig_id, .keep_all = TRUE)
}

read_one_sample <- function(sample) {
  contig_path <- file.path(
    options$megahit_dir,
    sample,
    "final.contigs.fa"
  )
  kraken_output_path <- file.path(options$kraken_dir, paste0(sample, ".out"))
  kraken_report_path <- file.path(options$kraken_dir, paste0(sample, ".report"))

  require_file(contig_path, paste0("MEGAHIT contigs for ", sample))
  require_file(kraken_output_path, paste0("Kraken2 output for ", sample))
  require_file(kraken_report_path, paste0("Kraken2 report for ", sample))

  sequences <- Biostrings::readDNAStringSet(contig_path)
  contig_lengths <- tibble::tibble(
    contig_id = sub("\\s.*$", "", names(sequences)),
    length_bp = as.numeric(Biostrings::width(sequences))
  )

  if (anyDuplicated(contig_lengths$contig_id)) {
    stop(sample, " FASTA contains duplicate contig IDs.", call. = FALSE)
  }

  kraken_output <- read_kraken_output(kraken_output_path)
  kraken_report <- read_kraken_report(kraken_report_path)

  joined <- contig_lengths |>
    dplyr::left_join(kraken_output, by = "contig_id") |>
    dplyr::left_join(kraken_report, by = "taxid") |>
    dplyr::mutate(
      classification = dplyr::coalesce(classification, "U"),
      taxid = dplyr::coalesce(taxid, "0"),
      taxon_name = dplyr::coalesce(
        taxon_name_inline,
        taxon_name_report,
        "Unclassified"
      ),
      rank = dplyr::coalesce(rank, "U"),
      sample = sample
    )

  missing_classifications <- sum(is.na(joined$taxon_name_inline) &
    is.na(joined$taxon_name_report))
  if (missing_classifications > 0L) {
    warning(
      sample,
      ": ",
      missing_classifications,
      " FASTA contig(s) had no resolved Kraken2 classification."
    )
  }

  length_disagreement <- joined |>
    dplyr::filter(
      !is.na(kraken_length),
      kraken_length != length_bp
    ) |>
    nrow()

  if (length_disagreement > 0L) {
    warning(
      sample,
      ": ",
      length_disagreement,
      " Kraken2 length(s) differ from FASTA lengths; FASTA lengths were retained."
    )
  }

  joined |>
    dplyr::select(
      sample,
      contig_id,
      length_bp,
      classification,
      taxid,
      taxon_name,
      rank
    )
}

all_contigs <- purrr::map_dfr(metadata$sample, read_one_sample) |>
  dplyr::left_join(metadata, by = "sample") |>
  dplyr::arrange(order, dplyr::desc(length_bp))

length_by_taxon <- all_contigs |>
  dplyr::group_by(
    sample,
    label,
    patient,
    timepoint,
    order,
    taxid,
    taxon_name,
    rank
  ) |>
  dplyr::summarise(
    contig_count = dplyr::n(),
    total_bp = sum(length_bp, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::mutate(total_kb = total_bp / 1000) |>
  dplyr::arrange(order, dplyr::desc(total_bp))

top_taxa <- length_by_taxon |>
  dplyr::group_by(sample) |>
  dplyr::slice_max(total_bp, n = top_n, with_ties = FALSE) |>
  dplyr::ungroup() |>
  dplyr::select(sample, taxon_name) |>
  dplyr::mutate(is_top = TRUE)

length_for_plot <- length_by_taxon |>
  dplyr::left_join(top_taxa, by = c("sample", "taxon_name")) |>
  dplyr::mutate(
    is_top = tidyr::replace_na(is_top, FALSE),
    taxon_display = dplyr::if_else(is_top, taxon_name, "Other")
  ) |>
  dplyr::group_by(sample, label, patient, timepoint, order, taxon_display) |>
  dplyr::summarise(
    contig_count = sum(contig_count),
    total_bp = sum(total_bp),
    total_kb = sum(total_kb),
    .groups = "drop"
  ) |>
  dplyr::arrange(order, dplyr::desc(total_bp))

readr::write_csv(all_contigs, file.path(table_dir, "contig_taxonomy.csv"))
readr::write_csv(
  length_by_taxon,
  file.path(table_dir, "length_by_taxon.csv")
)
readr::write_csv(
  length_for_plot,
  file.path(table_dir, "length_by_taxon_for_plot.csv")
)

ordered_labels <- metadata |>
  dplyr::arrange(order) |>
  dplyr::pull(label)

length_for_plot <- length_for_plot |>
  dplyr::mutate(label = factor(label, levels = ordered_labels))

taxon_levels <- c(
  sort(setdiff(unique(length_for_plot$taxon_display), "Other")),
  "Other"
)
non_other_taxa <- setdiff(taxon_levels, "Other")
taxon_colours <- if (length(non_other_taxa) > 0L) {
  stats::setNames(
    grDevices::hcl.colors(length(non_other_taxa), palette = "Dynamic"),
    non_other_taxa
  )
} else {
  character()
}
taxon_colours <- c(taxon_colours, Other = "grey70")

taxonomy_plot <- ggplot2::ggplot(
  length_for_plot,
  ggplot2::aes(x = label, y = total_kb, fill = taxon_display)
) +
  ggplot2::geom_col() +
  ggplot2::scale_fill_manual(values = taxon_colours, drop = FALSE) +
  ggplot2::labs(
    title = "Cumulative length of unmapped contigs by taxonomy",
    subtitle = paste0("Top ", top_n, " taxa selected separately for each sample"),
    x = "Sample",
    y = "Cumulative contig length (kb)",
    fill = "Taxonomy"
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
    panel.grid.major.x = ggplot2::element_blank(),
    legend.position = "right"
  )

ggplot2::ggsave(
  file.path(figure_dir, "unmapped_contig_taxonomy.png"),
  taxonomy_plot,
  width = 12,
  height = 7,
  dpi = 300
)
ggplot2::ggsave(
  file.path(figure_dir, "unmapped_contig_taxonomy.pdf"),
  taxonomy_plot,
  width = 12,
  height = 7
)

write_session_info(file.path(options$out_dir, "sessionInfo.txt"))

message(
  "Processed ",
  dplyr::n_distinct(all_contigs$sample),
  " samples. Outputs: ",
  normalizePath(options$out_dir, mustWork = TRUE)
)
