#' Fuse encoded modalities into meta-features
#'
#' @param mode one of 'early','late','hybrid'
#' @export
grf_fuse <- function(task, mode = c("early","late","hybrid")){
  stopifnot(inherits(task, "grf_task"))
  mode <- match.arg(mode)
  enc <- task$encodings
  if (length(enc) == 0) stop("call grf_encode() first")

  id_lists <- purrr::imap(enc, function(tbl, nm){
    if (!"..rowid.." %in% names(tbl)) {
      stop("Encoded modality '", nm, "' is missing the `..rowid..` column")
    }
    as.character(tbl$`..rowid..`)
  })

  common_ids <- id_lists[[1]]
  if (length(common_ids) == 0) {
    stop("No common sample identifiers across encoded modalities")
  }
  if (length(id_lists) > 1) {
    for (ids in id_lists[-1]) {
      common_ids <- common_ids[common_ids %in% ids]
    }
  }
  common_ids <- as.character(common_ids)
  if (length(common_ids) == 0) {
    stop("No common sample identifiers across encoded modalities")
  }

  aligned <- purrr::imap(enc, function(tbl, nm){
    grf_align_encoding(tbl, nm, common_ids)
  })

  combined <- NULL
  if (mode %in% c("early", "hybrid")) {
    prefixed <- purrr::imap(aligned, function(tbl, nm){
      feature_cols <- setdiff(names(tbl), "..rowid..")
      if (!length(feature_cols)) {
        return(tbl)
      }
      renamed <- tbl
      idx <- match(feature_cols, names(renamed))
      names(renamed)[idx] <- paste0(nm, "__", feature_cols)
      renamed
    })
    combined <- purrr::reduce(prefixed, function(acc, tbl){
      dplyr::full_join(acc, tbl, by = "..rowid..")
    })
    combined <- combined[, c("..rowid..", setdiff(names(combined), "..rowid..")), drop = FALSE]
  }

  task$fusion <- list(
    mode = mode,
    X = if (mode == "late") NULL else combined,
    per_modality = if (mode == "early") NULL else aligned
  )
  task$fit <- NULL
  task
}


grf_align_encoding <- function(tbl, modality, ids){
  if (!"..rowid.." %in% names(tbl)) {
    stop("Encoded modality '", modality, "' is missing the `..rowid..` column")
  }

  tbl <- tibble::as_tibble(tbl)
  tbl$`..rowid..` <- as.character(tbl$`..rowid..`)
  tbl <- tbl[tbl$`..rowid..` %in% ids, , drop = FALSE]

  ord <- match(ids, tbl$`..rowid..`)
  if (anyNA(ord)) {
    missing_ids <- ids[is.na(ord)]
    stop(
      "Encoded modality '", modality, "' is missing samples required for fusion: ",
      paste(missing_ids, collapse = ", ")
    )
  }

  aligned <- tbl[ord, , drop = FALSE]
  aligned <- aligned[, c("..rowid..", setdiff(names(aligned), "..rowid..")), drop = FALSE]
  tibble::as_tibble(aligned)
}
