#' Explain a fitted multimodal model using DALEX
#' @export
chm_explain <- function(task){
  stopifnot(inherits(task, "chm_task"))
  if (is.null(task$fit)) stop("Model not fitted")
  if (task$fit$type != "single") stop("Only single-model explanation in MVP")
  wf <- task$fit$fit
  fit <- workflows::extract_fit_parsnip(wf)$fit
  rec <- workflows::extract_recipe(wf)
  df <- recipes::juice(rec)
  y <- df[[task$fit$outcome]]
  x <- dplyr::select(df, -all_of(task$fit$outcome))
  f <- function(newdata){
    preds <- stats::predict(fit, newdata, type = if (task$fit$family=="classification") "response" else "response")
    as.numeric(preds)
  }
  DALEX::explain(model = fit, data = x, y = y, predict_function = f, label = "chimera")
}
