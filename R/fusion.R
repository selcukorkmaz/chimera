#' Fuse encoded modalities into meta-features
#'
#' @param mode one of 'early','late','hybrid'
#' @export
grf_fuse <- function(task, mode = c("early","late","hybrid")){
  stopifnot(inherits(task, "grf_task"))
  mode <- match.arg(mode)
  enc <- task$encodings
  if (length(enc) == 0) stop("call grf_encode() first")

  task_ids <- rownames(MultiAssayExperiment::colData(task$mae))
  aligned <- purrr::imap(enc, function(tbl, nm){
    grf_align_encoding(tbl, nm, task_ids)
  })

  combined <- purrr::reduce(aligned, function(acc, tbl){
    dplyr::full_join(acc, tbl, by = "..rowid..")
  })
  combined <- combined[, c("..rowid..", setdiff(names(combined), "..rowid..")), drop = FALSE]

  if (mode == "early"){
    task$fusion <- list(mode = mode, X = combined, per_modality = NULL)
  } else {
    task$fusion <- list(mode = mode, X = combined, per_modality = aligned)
  }
  task$fit <- NULL
  task
}


grf_align_encoding <- function(tbl, modality, ids){
  if (!"..rowid.." %in% names(tbl)) {
    stop("Encoded modality '", modality, "' is missing the `..rowid..` column")
  }
  tbl <- dplyr::mutate(tbl, `..rowid..` = as.character(`..rowid..`))
  feature_cols <- setdiff(names(tbl), "..rowid..")
  if (length(feature_cols)) {
    tbl <- dplyr::rename_with(tbl, ~paste0(modality, "__", .x), .cols = feature_cols)
  }
  scaffold <- tibble::tibble(`..rowid..` = ids)
  joined <- dplyr::left_join(scaffold, tbl, by = "..rowid..")
  tibble::as_tibble(joined)
}
