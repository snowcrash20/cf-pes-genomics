parse_cli_args <- function(args, defaults = list()) {
  if (length(args) == 0L) {
    return(defaults)
  }

  if (length(args) %% 2L != 0L) {
    stop(
      "Arguments must be supplied as '--name value' pairs. Use --help for usage.",
      call. = FALSE
    )
  }

  keys <- args[seq.int(1L, length(args), by = 2L)]
  values <- args[seq.int(2L, length(args), by = 2L)]

  if (any(!grepl("^--[a-z0-9-]+$", keys))) {
    stop(
      "Every argument name must begin with '--'. Use --help for usage.",
      call. = FALSE
    )
  }

  keys <- sub("^--", "", keys)
  keys <- gsub("-", "_", keys, fixed = TRUE)

  for (index in seq_along(keys)) {
    defaults[[keys[[index]]]] <- values[[index]]
  }

  defaults
}

require_cli_values <- function(options, names) {
  missing <- names[vapply(
    names,
    function(name) {
      is.null(options[[name]]) ||
        length(options[[name]]) == 0L ||
        is.na(options[[name]]) ||
        !nzchar(options[[name]])
    },
    logical(1)
  )]

  if (length(missing) > 0L) {
    stop(
      "Missing required argument(s): ",
      paste(paste0("--", gsub("_", "-", missing, fixed = TRUE)), collapse = ", "),
      ". Use --help for usage.",
      call. = FALSE
    )
  }
}

require_packages <- function(packages) {
  unavailable <- packages[!vapply(
    packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )]

  if (length(unavailable) > 0L) {
    stop(
      "Missing R package(s): ",
      paste(unavailable, collapse = ", "),
      ". Run Rscript scripts/install_r_dependencies.R, then retry.",
      call. = FALSE
    )
  }
}

require_file <- function(path, description = "File") {
  if (!file.exists(path)) {
    stop(description, " not found: ", path, call. = FALSE)
  }
}

require_directory <- function(path, description = "Directory") {
  if (!dir.exists(path)) {
    stop(description, " not found: ", path, call. = FALSE)
  }
}

parse_number <- function(value, name, minimum = -Inf, maximum = Inf) {
  parsed <- suppressWarnings(as.numeric(value))

  if (
    length(parsed) != 1L ||
      is.na(parsed) ||
      parsed < minimum ||
      parsed > maximum
  ) {
    stop(
      "--",
      gsub("_", "-", name, fixed = TRUE),
      " must be a number between ",
      minimum,
      " and ",
      maximum,
      ".",
      call. = FALSE
    )
  }

  parsed
}

parse_integer <- function(value, name, minimum = 1L) {
  parsed <- suppressWarnings(as.integer(value))

  if (length(parsed) != 1L || is.na(parsed) || parsed < minimum) {
    stop(
      "--",
      gsub("_", "-", name, fixed = TRUE),
      " must be an integer greater than or equal to ",
      minimum,
      ".",
      call. = FALSE
    )
  }

  parsed
}

read_sample_metadata <- function(path) {
  require_file(path, "Sample metadata")

  metadata <- readr::read_tsv(
    path,
    show_col_types = FALSE,
    progress = FALSE,
    col_types = readr::cols(.default = readr::col_character())
  )

  required_columns <- c("sample", "label", "patient", "timepoint", "order")
  missing_columns <- setdiff(required_columns, names(metadata))

  if (length(missing_columns) > 0L) {
    stop(
      "Sample metadata is missing column(s): ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  if (anyDuplicated(metadata$sample)) {
    duplicate_samples <- unique(metadata$sample[duplicated(metadata$sample)])
    stop(
      "Sample metadata contains duplicate sample IDs: ",
      paste(duplicate_samples, collapse = ", "),
      call. = FALSE
    )
  }

  metadata <- dplyr::mutate(
    metadata,
    order = suppressWarnings(as.integer(order))
  )

  if (any(is.na(metadata$order))) {
    stop("Every metadata 'order' value must be an integer.", call. = FALSE)
  }

  metadata
}

slugify <- function(value) {
  value <- tolower(value)
  value <- gsub("[^a-z0-9]+", "-", value)
  value <- gsub("(^-|-$)", "", value)
  value
}

write_session_info <- function(path) {
  lines <- capture.output(utils::sessionInfo())
  writeLines(lines, con = path, useBytes = TRUE)
}
