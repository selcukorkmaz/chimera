#' Create a multimodal task
#'
#' @param mae MultiAssayExperiment or list of modality tables. For MVP we accept
#' either; lists will be coerced to MAE with synthetic assays.
#' @param outcome A tibble/data.frame with one column outcome (factor/numeric) and optional id column name `sample_id`.
#' @param id_col Name of sample id column shared across modalities (default 'sample_id').
#' @return An object of class `grf_task`.
#' @export
grf_task <- function(mae, outcome, id_col = "sample_id"){
  outcome_df <- as.data.frame(outcome, stringsAsFactors = FALSE)
  if (nrow(outcome_df) == 0) {
    stop("Outcome must contain at least one sample")
  }

  if (id_col %in% names(outcome_df)) {
    sample_ids <- outcome_df[[id_col]]
  } else if (!is.null(rownames(outcome_df))) {
    sample_ids <- rownames(outcome_df)
    outcome_df[[id_col]] <- sample_ids
  } else {
    stop("Outcome must include a sample identifier column or rownames")
  }

  sample_ids <- as.character(sample_ids)
  if (anyNA(sample_ids)) stop("Sample identifiers cannot contain NA")
  if (anyDuplicated(sample_ids)) stop("Sample identifiers must be unique")
  rownames(outcome_df) <- sample_ids

  normalize_se <- function(se) {
    assays_list <- SummarizedExperiment::assays(se)
    if (length(assays_list) == 0) {
      stop("Each modality must contain at least one assay")
    }
    normalized <- lapply(assays_list, function(mat) {
      mat <- as.matrix(mat)
      matches_col <- if (!is.null(colnames(mat))) sum(colnames(mat) %in% sample_ids) else 0
      matches_row <- if (!is.null(rownames(mat))) sum(rownames(mat) %in% sample_ids) else 0
      if ((matches_row > matches_col) || (matches_col == 0 && !is.null(rownames(mat)) && is.null(colnames(mat)))) {
        mat <- t(mat)
      }
      if (is.null(colnames(mat))) {
        stop("Assays must have sample identifiers in either columns or rows")
      }
      mat
    })
    names(normalized) <- names(assays_list)
    SummarizedExperiment::SummarizedExperiment(assays = normalized)
  }

  se_from_matrix <- function(mat) {
    mat <- as.matrix(mat)
    matches_col <- if (!is.null(colnames(mat))) sum(colnames(mat) %in% sample_ids) else 0
    matches_row <- if (!is.null(rownames(mat))) sum(rownames(mat) %in% sample_ids) else 0
    if ((matches_row > matches_col) || (matches_col == 0 && !is.null(rownames(mat)) && is.null(colnames(mat)))) {
      mat <- t(mat)
    }
    if (is.null(colnames(mat))) {
      stop("Matrix modalities require sample identifiers in dimnames")
    }
    SummarizedExperiment::SummarizedExperiment(assays = list(data = mat))
  }

  se_from_df <- function(df) {
    df <- as.data.frame(df, stringsAsFactors = FALSE)
    if (id_col %in% names(df)) {
      ids <- df[[id_col]]
      df[[id_col]] <- NULL
    } else if (!is.null(rownames(df))) {
      ids <- rownames(df)
    } else {
      stop("Each modality table must include the sample id column or rownames")
    }
    ids <- as.character(ids)
    if (anyNA(ids)) stop("Sample identifiers cannot contain NA")
    if (anyDuplicated(ids)) stop("Sample identifiers must be unique")
    mat <- t(as.matrix(df))
    colnames(mat) <- ids
    SummarizedExperiment::SummarizedExperiment(assays = list(data = mat))
  }

  if (inherits(mae, "MultiAssayExperiment")) {
    experiments_raw <- as.list(MultiAssayExperiment::experiments(mae))
    experiments <- lapply(experiments_raw, normalize_se)
  } else if (is.list(mae)) {
    if (length(mae) == 0) stop("No modalities supplied")
    if (is.null(names(mae))) names(mae) <- paste0("modality", seq_along(mae))
    experiments <- lapply(mae, function(x) {
      if (inherits(x, "SummarizedExperiment")) {
        normalize_se(x)
      } else if (is.matrix(x) || inherits(x, "Matrix")) {
        se_from_matrix(x)
      } else if (is.data.frame(x)) {
        se_from_df(x)
      } else {
        stop("Each modality must be a SummarizedExperiment, matrix, or data.frame")
      }
    })
  } else stop("mae must be MultiAssayExperiment or list of modalities")

  if (is.null(names(experiments))) {
    names(experiments) <- paste0("modality", seq_along(experiments))
  }

  experiment_ids <- lapply(experiments, function(se) {
    ids <- colnames(SummarizedExperiment::assay(se))
    if (is.null(ids)) stop("Assays must contain column names for samples")
    ids
  })

  common_ids <- sample_ids
  for (ids in experiment_ids) {
    common_ids <- common_ids[common_ids %in% ids]
  }
  common_ids <- unique(common_ids)
  if (length(common_ids) == 0) {
    stop("No overlapping sample identifiers between outcome and modalities")
  }

  outcome_aligned <- outcome_df[common_ids, , drop = FALSE]

  experiments <- lapply(experiments, function(se) {
    assays_list <- SummarizedExperiment::assays(se)
    subsetted <- lapply(assays_list, function(mat) {
      mat[, common_ids, drop = FALSE]
    })
    names(subsetted) <- names(assays_list)
    SummarizedExperiment::SummarizedExperiment(assays = subsetted)
  })

  mae_obj <- MultiAssayExperiment::MultiAssayExperiment(
    experiments = MultiAssayExperiment::ExperimentList(experiments),
    colData = S4Vectors::DataFrame(outcome_aligned)
  )

  structure(list(
    mae = mae_obj,
    id_col = id_col,
    modalities = names(MultiAssayExperiment::experiments(mae_obj)),
    registry = list(),
    encodings = list(),
    fusion = NULL,
    fit = NULL
  ), class = "grf_task")
}

#' @export
print.grf_task <- function(x, ...){
  stopifnot(inherits(x, "grf_task"))
  cd <- MultiAssayExperiment::colData(x$mae)
  n_samples <- nrow(cd)
  cat("<grf_task>\n")
  cat("  Samples: ", n_samples, "\n", sep = "")
  if (length(x$modalities)) {
    cat("  Modalities (", length(x$modalities), "): ", paste(x$modalities, collapse = ", "), "\n", sep = "")
  } else {
    cat("  Modalities: none\n")
  }
  outcome_cols <- setdiff(colnames(cd), x$id_col)
  if (length(outcome_cols)) {
    cat("  Outcome columns: ", paste(outcome_cols, collapse = ", "), "\n", sep = "")
  } else {
    cat("  Outcome columns: <none>\n")
  }
  invisible(x)
}


#' Register a modality encoder
#'
#' @param task grf_task
#' @param name modality name present in task$modalities
#' @param encoder a function taking (se, id_col) and returning a feature matrix/data.frame with rownames as sample ids
#' @export
grf_add_modality <- function(task, name, encoder = c("numeric","tabular","text")) {
  stopifnot(inherits(task, "grf_task"))
  if (!name %in% task$modalities) stop("Unknown modality: ", name)

  if (is.character(encoder)) {
    encoder <- match.arg(encoder)
    enc_fun <- switch(
      encoder,
      numeric = grf_encoder_numeric(),
      tabular = grf_encoder_tabular(),
      text = grf_encoder_text(),
      stop("Unknown encoder keyword")
    )
  } else if (is.function(encoder)) {
    enc_fun <- encoder
  } else {
    stop("encoder must be a keyword or function")
  }

  task$registry[[name]] <- enc_fun
  task
}


#' Encode all registered modalities
#' @export
grf_encode <- function(task){
  stopifnot(inherits(task, "grf_task"))
  if (length(task$registry) == 0) stop("No modality encoders registered; use grf_add_modality() first")

  task_ids <- rownames(MultiAssayExperiment::colData(task$mae))
  enc <- purrr::imap(task$registry, function(enc_fun, nm){
    se <- MultiAssayExperiment::experiments(task$mae)[[nm]]
    if (is.null(se)) {
      stop("Modality not found in task: ", nm)
    }
    encoded <- enc_fun(se, task$id_col)
    encoded_tbl <- grf_validate_encoding_output(encoded, nm)
    missing_ids <- setdiff(encoded_tbl$`..rowid..`, task_ids)
    if (length(missing_ids)) {
      stop(
        "Encoder for modality '", nm, "' returned sample ids not present in the task: ",
        paste(missing_ids, collapse = ", ")
      )
    }
    encoded_tbl
  })
  task$encodings <- enc
  task$fusion <- NULL
  task$fit <- NULL
  task
}


# Internal helpers ---------------------------------------------------------

grf_encoder_numeric <- function(pca_k = NULL){
  if (!is.null(pca_k)) {
    if (!is.numeric(pca_k) || length(pca_k) != 1 || is.na(pca_k) || pca_k <= 0) {
      stop("`pca_k` must be a positive integer when supplied")
    }
    pca_k <- as.integer(pca_k)
  }

  function(se, id_col){
    samples <- grf_modality_samples(se)
    if (!is.numeric(samples)) {
      stop("Numeric encoder requires numeric assay values")
    }

    sample_ids <- rownames(samples)
    if (is.null(sample_ids)) {
      stop("Numeric modality is missing sample identifiers")
    }

    if (ncol(samples) == 0) {
      return(grf_matrix_to_tibble(grf_empty_feature_matrix(sample_ids)))
    }

    scaled <- grf_scale_matrix(samples)

    if (!is.null(pca_k)) {
      max_rank <- min(ncol(scaled), nrow(scaled) - 1)
      if (!is.na(max_rank) && max_rank > 0) {
        k <- min(pca_k, max_rank)
        scaled_for_pca <- scaled
        scaled_for_pca[is.na(scaled_for_pca)] <- 0
        pca <- stats::prcomp(scaled_for_pca, center = FALSE, scale. = FALSE, rank. = k)
        pcs <- pca$x
        if (!is.null(pcs) && ncol(pcs) > 0) {
          k <- min(k, ncol(pcs))
          pcs <- pcs[, seq_len(k), drop = FALSE]
          colnames(pcs) <- paste0("PC", seq_len(ncol(pcs)))
          return(grf_matrix_to_tibble(pcs))
        }
      }
    }

    grf_matrix_to_tibble(scaled)
  }
}


grf_encoder_tabular <- function(){
  function(se, id_col){
    samples <- grf_modality_samples(se)
    sample_ids <- rownames(samples)
    if (is.null(sample_ids)) {
      stop("Tabular modality is missing sample identifiers")
    }

    df <- as.data.frame(samples, stringsAsFactors = FALSE)
    if (!ncol(df)) {
      return(grf_matrix_to_tibble(grf_empty_feature_matrix(sample_ids)))
    }

    numeric_cols <- list()
    cat_mats <- list()

    for (nm in names(df)) {
      prepared <- grf_tabular_prepare_column(df[[nm]])
      if (prepared$type == "numeric") {
        numeric_cols[[nm]] <- prepared$values
      } else if (prepared$type == "categorical") {
        dummy <- grf_tabular_one_hot(prepared$values, sample_ids, nm)
        if (!is.null(dummy)) {
          cat_mats[[length(cat_mats) + 1L]] <- dummy
        }
      }
    }

    matrices <- list()
    if (length(numeric_cols)) {
      numeric_df <- as.data.frame(numeric_cols, stringsAsFactors = FALSE)
      rownames(numeric_df) <- sample_ids
      numeric_mat <- as.matrix(numeric_df)
      if (ncol(numeric_mat) > 0) {
        matrices <- c(matrices, list(grf_scale_matrix(numeric_mat)))
      }
    }
    if (length(cat_mats)) {
      matrices <- c(matrices, cat_mats)
    }

    if (!length(matrices)) {
      return(grf_matrix_to_tibble(grf_empty_feature_matrix(sample_ids)))
    }

    combined <- do.call(cbind, matrices)
    rownames(combined) <- sample_ids
    grf_matrix_to_tibble(combined)
  }
}


grf_encoder_text <- function(max_features = NULL, term_count_min = NULL, doc_prop_max = NULL){
  if (!is.null(max_features)) {
    if (!is.numeric(max_features) || length(max_features) != 1 || is.na(max_features) || max_features <= 0) {
      stop("`max_features` must be a positive integer when supplied")
    }
    max_features <- as.integer(max_features)
  }
  if (!is.null(term_count_min)) {
    if (!is.numeric(term_count_min) || length(term_count_min) != 1 || is.na(term_count_min) || term_count_min <= 0) {
      stop("`term_count_min` must be a positive integer when supplied")
    }
    term_count_min <- as.integer(term_count_min)
  }
  if (!is.null(doc_prop_max)) {
    if (!is.numeric(doc_prop_max) || length(doc_prop_max) != 1 || is.na(doc_prop_max) || doc_prop_max <= 0 || doc_prop_max > 1) {
      stop("`doc_prop_max` must be a numeric value in (0, 1]")
    }
    doc_prop_max <- as.numeric(doc_prop_max)
  }

  function(se, id_col){
    samples <- grf_modality_samples(se)
    sample_ids <- rownames(samples)
    if (is.null(sample_ids)) stop("Text encoder requires sample identifiers")

    df <- as.data.frame(samples, stringsAsFactors = FALSE)
    text_cols <- names(df)[tolower(names(df)) == "text"]
    if (!length(text_cols)) {
      stop("Text modality must contain at least one column named 'text'")
    }

    text_df <- df[text_cols, drop = FALSE]
    text_vec <- vapply(seq_len(nrow(text_df)), function(i){
      row <- text_df[i, , drop = TRUE]
      row <- row[!is.na(row)]
      row <- trimws(as.character(row))
      row <- row[nzchar(row)]
      if (length(row)) paste(row, collapse = " ") else ""
    }, character(1))
    names(text_vec) <- sample_ids

    tokens <- text2vec::itoken(text_vec, ids = sample_ids, progressbar = FALSE)
    vocab <- text2vec::create_vocabulary(tokens)

    prune_args <- list(vocab = vocab)
    if (!is.null(term_count_min)) prune_args$term_count_min <- term_count_min
    if (!is.null(doc_prop_max)) prune_args$doc_proportion_max <- doc_prop_max
    if (!is.null(max_features)) prune_args$max_number_of_terms <- max_features
    if (length(prune_args) > 1) {
      vocab <- do.call(text2vec::prune_vocabulary, prune_args)
    }

    if (nrow(vocab) == 0) {
      return(grf_matrix_to_tibble(grf_empty_feature_matrix(sample_ids)))
    }

    vectorizer <- text2vec::vocab_vectorizer(vocab)
    tokens <- text2vec::itoken(text_vec, ids = sample_ids, progressbar = FALSE)
    dtm <- text2vec::create_dtm(tokens, vectorizer)
    tfidf <- text2vec::TfIdf$new()
    tfidf_mat <- tfidf$fit_transform(dtm)

    if (ncol(tfidf_mat) == 0) {
      return(grf_matrix_to_tibble(grf_empty_feature_matrix(sample_ids)))
    }

    dense <- as.matrix(tfidf_mat)
    rownames(dense) <- sample_ids
    if (ncol(dense) > 0) {
      colnames(dense) <- make.unique(paste0("text_", colnames(dense)))
    }
    grf_matrix_to_tibble(dense)
  }
}


grf_modality_samples <- function(se){
  mat <- SummarizedExperiment::assay(se)
  if (is.null(mat)) stop("Modality assay is empty")
  mat <- as.matrix(mat)
  sample_ids <- colnames(mat)
  if (is.null(sample_ids)) stop("Assay must provide sample identifiers in column names")
  mat <- t(mat)
  rownames(mat) <- sample_ids
  mat
}


grf_matrix_to_tibble <- function(mat){
  if (inherits(mat, "Matrix")) {
    mat <- as.matrix(mat)
  }
  if (!is.matrix(mat)) {
    stop("Encoder output must be matrix-like")
  }
  if (ncol(mat) > 0) {
    if (is.null(colnames(mat))) {
      colnames(mat) <- paste0("feature_", seq_len(ncol(mat)))
    } else {
      colnames(mat) <- make.unique(colnames(mat))
    }
  }
  if (is.null(rownames(mat))) {
    stop("Matrix must contain rownames for sample identifiers")
  }
  df <- tibble::rownames_to_column(as.data.frame(mat, stringsAsFactors = FALSE), "..rowid..")
  tibble::as_tibble(df)
}


grf_empty_feature_matrix <- function(sample_ids){
  sample_ids <- as.character(sample_ids)
  empty <- matrix(nrow = length(sample_ids), ncol = 0)
  rownames(empty) <- sample_ids
  empty
}


grf_scale_matrix <- function(mat){
  if (inherits(mat, "Matrix")) {
    mat <- as.matrix(mat)
  }
  if (!is.matrix(mat)) {
    stop("Expected a matrix for scaling")
  }
  rn <- rownames(mat)
  cn <- colnames(mat)
  storage.mode(mat) <- "double"
  if (!ncol(mat)) {
    rownames(mat) <- rn
    colnames(mat) <- cn
    return(mat)
  }
  centers <- apply(mat, 2, function(x){
    if (all(is.na(x))) return(0)
    mean(x, na.rm = TRUE)
  })
  centers[is.nan(centers)] <- 0
  scaled <- sweep(mat, 2, centers, "-")
  scales <- apply(mat, 2, stats::sd, na.rm = TRUE)
  scales[is.na(scales) | scales == 0] <- 1
  scaled <- sweep(scaled, 2, scales, "/")
  scaled[is.na(mat)] <- NA_real_
  rownames(scaled) <- rn
  colnames(scaled) <- cn
  scaled
}


grf_tabular_prepare_column <- function(x){
  if (inherits(x, c("Date", "POSIXct", "POSIXlt", "POSIXt"))) {
    numeric_vals <- as.numeric(x)
    numeric_vals[is.na(x)] <- NA_real_
    return(list(type = "numeric", values = numeric_vals))
  }
  if (is.factor(x)) {
    x <- as.character(x)
  }
  if (is.logical(x)) {
    numeric_vals <- as.numeric(x)
    numeric_vals[is.na(x)] <- NA_real_
    return(list(type = "numeric", values = numeric_vals))
  }
  if (is.numeric(x)) {
    return(list(type = "numeric", values = as.numeric(x)))
  }
  if (is.character(x)) {
    trimmed <- trimws(x)
    trimmed[trimmed == ""] <- NA_character_
    suppressWarnings(num <- as.numeric(trimmed))
    non_missing <- !is.na(trimmed)
    if (any(non_missing) && all(!is.na(num[non_missing]))) {
      num[!non_missing] <- NA_real_
      return(list(type = "numeric", values = num))
    }
    if (!any(non_missing)) {
      return(list(type = "numeric", values = rep(NA_real_, length(trimmed))))
    }
    return(list(type = "categorical", values = trimmed))
  }
  char_vals <- as.character(x)
  char_vals <- trimws(char_vals)
  char_vals[char_vals == ""] <- NA_character_
  list(type = "categorical", values = char_vals)
}


grf_tabular_one_hot <- function(values, sample_ids, prefix){
  values <- as.character(values)
  values <- trimws(values)
  values[values == ""] <- NA_character_
  encoded <- values
  if (anyNA(encoded)) {
    encoded[is.na(encoded)] <- "__missing__"
  }
  levels <- unique(encoded)
  if (!length(levels)) {
    return(NULL)
  }
  mat <- matrix(0, nrow = length(encoded), ncol = length(levels))
  for (i in seq_along(levels)) {
    mat[, i] <- as.numeric(encoded == levels[i])
  }
  clean_levels <- levels
  clean_levels[clean_levels == "__missing__"] <- "missing"
  colnames(mat) <- make.unique(paste0(prefix, "__", make.names(clean_levels, unique = FALSE)))
  rownames(mat) <- sample_ids
  mat
}


grf_validate_encoding_output <- function(enc, modality){
  if (inherits(enc, "Matrix") || is.matrix(enc)) {
    enc_tbl <- grf_matrix_to_tibble(enc)
  } else if (tibble::is_tibble(enc)) {
    enc_tbl <- enc
  } else if (is.data.frame(enc)) {
    enc_tbl <- tibble::as_tibble(enc)
  } else {
    stop("Encoder for modality '", modality, "' must return a tibble, data.frame, or matrix")
  }

  if (!"..rowid.." %in% names(enc_tbl)) {
    stop("Encoder for modality '", modality, "' must include a `..rowid..` column")
  }

  enc_tbl <- dplyr::mutate(enc_tbl, `..rowid..` = as.character(`..rowid..`))
  if (anyNA(enc_tbl$`..rowid..`)) {
    stop("Encoder for modality '", modality, "' produced missing sample identifiers")
  }
  if (anyDuplicated(enc_tbl$`..rowid..`)) {
    stop("Encoder for modality '", modality, "' must return unique sample identifiers")
  }

  feature_cols <- setdiff(names(enc_tbl), "..rowid..")
  if (anyDuplicated(feature_cols)) {
    stop("Encoder for modality '", modality, "' returned duplicated feature names")
  }
  enc_tbl <- enc_tbl[, c("..rowid..", feature_cols), drop = FALSE]
  tibble::as_tibble(enc_tbl)
}
