#' Create a multimodal task
#'
#' @param mae MultiAssayExperiment or list of modality tables. For MVP we accept
#' either; lists will be coerced to MAE with synthetic assays.
#' @param outcome A tibble/data.frame with one column outcome (factor/numeric) and optional id column name `sample_id`.
#' @param id_col Name of sample id column shared across modalities (default 'sample_id').
#' @return An object of class `grf_task`.
#' @export
grf_task <- function(mae, outcome, id_col = "sample_id"){
  if (inherits(mae, "MultiAssayExperiment")) {
    mae_obj <- mae
  } else if (is.list(mae)) {
    # Coerce simple list of data.frames/matrices to MAE
    assays <- lapply(mae, function(x){
      if (is.data.frame(x)) x <- as.matrix(x)
      SummarizedExperiment::SummarizedExperiment(assays = list(counts = x))
    })
    mae_obj <- MultiAssayExperiment::MultiAssayExperiment(
      experiments = MultiAssayExperiment::ExperimentList(assays),
      colData = S4Vectors::DataFrame(outcome)
    )
  } else stop("mae must be MultiAssayExperiment or list of modalities")


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


#' Register a modality encoder
#'
#' @param task grf_task
#' @param name modality name present in task$modalities
#' @param encoder a function taking (se, id_col) and returning a feature matrix/data.frame with rownames as sample ids
#' @export
grf_add_modality <- function(task, name, encoder = c("numeric","text","image")) {
  stopifnot(inherits(task, "grf_task"))
  if (!name %in% task$modalities) stop("Unknown modality: ", name)

  if (is.character(encoder)) {
    enc_fun <- switch(encoder,
                      numeric = function(se, id_col) {
                        X <- SummarizedExperiment::assay(se)
                        X <- scale(as.matrix(X))
                        tibble::as_tibble(as.data.frame(X), rownames = "..rowid..")
                      },
                      text = function(se, id_col) {
                        # tf-idf code here
                      },
                      image = function(se, id_col) {
                        # image stats here
                      },
                      stop("Unknown encoder type")
    )
  } else {
    enc_fun <- encoder
  }

  task$registry[[name]] <- enc_fun
  task
}


#' Encode all registered modalities
#' @export
grf_encode <- function(task){
  stopifnot(inherits(task, "grf_task"))
  enc <- purrr::imap(task$registry, function(enc_fun, nm){
    se <- MultiAssayExperiment::experiments(task$mae)[[nm]]
    enc_fun(se, task$id_col)
  })
  task$encodings <- enc
  task
}

