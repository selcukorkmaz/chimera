#' Fuse encoded modalities into meta-features
#'
#' @param mode one of 'early','late','hybrid'
#' @export
grf_fuse <- function(task, mode = c("early","late","hybrid")){
  stopifnot(inherits(task, "grf_task"))
  mode <- match.arg(mode)
  enc <- task$encodings
  if (length(enc) == 0) stop("call grf_encode() first")


  # Align by id
  ids <- Reduce(intersect, lapply(enc, function(x) x$`..rowid..`))
  enc <- lapply(enc, function(x) dplyr::filter(x, `..rowid..` %in% ids))


  if (mode == "early"){
    X <- Reduce(function(a,b) dplyr::full_join(a,b, by = "..rowid.."), enc)
    task$fusion <- list(mode = mode, X = X, stacks = NULL)
  } else if (mode == "late"){
    task$fusion <- list(mode = mode, per_modality = enc, meta = NULL)
  } else if (mode == "hybrid"){
    X <- Reduce(function(a,b) dplyr::full_join(a,b, by = "..rowid.."), enc)
    task$fusion <- list(mode = mode, X = X, per_modality = enc)
  }
  task
}
