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
      numeric = grf_encoder_numeric,
      tabular = grf_encoder_tabular,
      text = grf_encoder_text,
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

grf_encoder_numeric <- function(se, id_col){
  samples <- grf_modality_samples(se)
  if (!is.numeric(samples)) {
    stop("Numeric encoder requires numeric assay values")
  }
  scaled <- scale(samples)
  scaled[is.na(scaled)] <- 0
  grf_matrix_to_tibble(scaled)
}


grf_encoder_tabular <- function(se, id_col){
  samples <- grf_modality_samples(se)
  grf_matrix_to_tibble(samples)
}


grf_encoder_text <- function(se, id_col){
  samples <- grf_modality_samples(se)
  sample_ids <- rownames(samples)
  if (is.null(sample_ids)) stop("Text encoder requires sample identifiers")

  text_vec <- vapply(seq_along(sample_ids), function(i){
    row <- samples[i, , drop = TRUE]
    row <- row[!is.na(row)]
    row <- trimws(as.character(row))
    row <- row[nzchar(row)]
    if (length(row)) paste(row, collapse = " ") else ""
  }, character(1), USE.NAMES = FALSE)
  names(text_vec) <- sample_ids

  tokens <- text2vec::itoken(text_vec, ids = sample_ids, progressbar = FALSE)
  vocab <- text2vec::create_vocabulary(tokens)
  if (nrow(vocab) == 0) {
    return(tibble::tibble(`..rowid..` = sample_ids))
  }
  vectorizer <- text2vec::vocab_vectorizer(vocab)
  tokens <- text2vec::itoken(text_vec, ids = sample_ids, progressbar = FALSE)
  dtm <- text2vec::create_dtm(tokens, vectorizer)
  tfidf <- text2vec::TfIdf$new()
  tfidf_mat <- tfidf$fit_transform(dtm)
  if (ncol(tfidf_mat) == 0) {
    return(tibble::tibble(`..rowid..` = sample_ids))
  }
  features <- grf_matrix_to_tibble(tfidf_mat)
  if (ncol(features) > 1) {
    feature_cols <- names(features)[-1]
    names(features)[-1] <- paste0("text_", feature_cols)
  }
  features
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
