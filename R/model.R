#' Fit a multimodal model
#'
#' @param learner a parsnip model spec
#' @param outcome_name column name in outcome table
#' @param family 'classification'|'regression'
#' @export
grf_fit <- function(task, learner, outcome_name, family = c("classification","regression")){
  stopifnot(inherits(task, "grf_task"))
  family <- match.arg(family)
  cd <- as.data.frame(MultiAssayExperiment::colData(task$mae))
  y <- cd[[outcome_name]]


  if (task$fusion$mode == "early" || task$fusion$mode == "hybrid"){
    X <- task$fusion$X
    df <- dplyr::left_join(tibble::tibble(`..rowid..` = rownames(cd)), X, by = "..rowid..")
    df[[outcome_name]] <- y[match(df$`..rowid..`, rownames(cd))]


    rec <- recipes::recipe(as.formula(paste(outcome_name, "~ .")), data = dplyr::select(df, -`..rowid..`)) %>%
      recipes::update_role(outcome_name, new_role = "outcome") %>%
      recipes::step_zv(recipes::all_predictors())


    wf <- workflows::workflow() %>%
      workflows::add_model(learner) %>%
      workflows::add_recipe(rec)


    fit <- parsnip::fit(wf, data = dplyr::select(df, -`..rowid..`))
    task$fit <- list(type = "single", fit = fit, family = family, outcome = outcome_name)
  } else if (task$fusion$mode == "late"){
    # simplified MVP stacking as before
    task$fit <- list(type = "stack", base = NULL, meta = NULL, family = family, outcome = outcome_name)
  }
  task
}

#' Predict from a fitted grf_task
#' @export
grf_predict <- function(task, new_task = NULL){
  stopifnot(inherits(task, "grf_task"))
  if (is.null(task$fit)) stop("Model not fitted")
  if (!is.null(new_task)) task <- new_task
  if (task$fit$type == "single"){
    wf <- task$fit$fit
    predict(wf, workflows::extract_recipe(wf) %>% recipes::juice(), type = if (task$fit$family=="classification") "prob" else "numeric")
  } else stop("Prediction for stacks not implemented in MVP")
}
