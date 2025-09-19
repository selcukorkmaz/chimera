local_omnibus_task <- function(modality_list, outcome){
  grf_task(modality_list, outcome)
}

test_that("numeric encoder scales features and supports PCA", {
  set.seed(123)
  sample_ids <- paste0("S", seq_len(40))
  feature_ids <- paste0("g", seq_len(15))
  mat <- matrix(rnorm(length(sample_ids) * length(feature_ids)),
                nrow = length(feature_ids), ncol = length(sample_ids),
                dimnames = list(feature_ids, sample_ids))
  se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = mat))
  modalities <- list(omics = se)
  outcome <- data.frame(sample_id = sample_ids,
                        status = sample(c("A", "B"), length(sample_ids), replace = TRUE))

  task <- local_omnibus_task(modalities, outcome)
  task <- grf_add_modality(task, "omics", "numeric")
  task <- grf_encode(task)
  encoded <- task$encodings$omics

  expect_s3_class(encoded, "tbl_df")
  expect_equal(nrow(encoded), length(sample_ids))
  expect_equal(ncol(encoded), length(feature_ids) + 1)

  feature_mat <- as.matrix(encoded[, -1, drop = FALSE])
  expect_true(all(abs(colMeans(feature_mat)) < 1e-6))
  expect_true(all(abs(apply(feature_mat, 2, stats::sd) - 1) < 1e-6))

  task_pca <- local_omnibus_task(modalities, outcome)
  task_pca <- grf_add_modality(task_pca, "omics", grf_encoder_numeric(pca_k = 5))
  task_pca <- grf_encode(task_pca)
  encoded_pca <- task_pca$encodings$omics

  expect_equal(ncol(encoded_pca), 6)
  expect_equal(names(encoded_pca)[-1], paste0("PC", seq_len(5)))
})

test_that("tabular encoder scales numeric columns and one-hot encodes categoricals", {
  sample_ids <- paste0("P", seq_len(8))
  clinical <- data.frame(
    sample_id = sample_ids,
    age = seq(45, 59, length.out = length(sample_ids)),
    bmi = as.character(round(seq(21, 28, length.out = length(sample_ids)), 1)),
    sex = rep(c("F", "M"), length.out = length(sample_ids)),
    smoker = factor(rep(c("yes", "no"), length.out = length(sample_ids))),
    stringsAsFactors = FALSE
  )
  outcome <- data.frame(sample_id = sample_ids, status = rnorm(length(sample_ids)))

  task <- local_omnibus_task(list(clinical = clinical), outcome)
  task <- grf_add_modality(task, "clinical", "tabular")
  task <- grf_encode(task)

  encoded <- task$encodings$clinical
  expect_s3_class(encoded, "tbl_df")
  expect_equal(nrow(encoded), length(sample_ids))

  expect_true("age" %in% names(encoded))
  expect_true("bmi" %in% names(encoded))
  expect_true(all(abs(mean(encoded$age)) < 1e-8))
  expect_true(all(abs(stats::sd(encoded$age) - 1) < 1e-8))
  expect_true(all(abs(mean(encoded$bmi)) < 1e-8))
  expect_true(all(abs(stats::sd(encoded$bmi) - 1) < 1e-8))

  sex_cols <- grep("^sex__", names(encoded), value = TRUE)
  smoker_cols <- grep("^smoker__", names(encoded), value = TRUE)
  expect_true(length(sex_cols) >= 1)
  expect_true(length(smoker_cols) >= 1)
  for (col in c(sex_cols, smoker_cols)) {
    expect_true(all(encoded[[col]] %in% c(0, 1)))
  }
})

test_that("text encoder creates tf-idf features with optional limits", {
  sample_ids <- paste0("S", seq_len(6))
  text_tbl <- data.frame(
    sample_id = sample_ids,
    text = c(
      "cancer growth signal",
      "normal tissue baseline",
      "tumor cell growth",
      "benign tissue sample",
      "cancer cell mutation",
      "normal sample analysis"
    ),
    stringsAsFactors = FALSE
  )
  outcome <- data.frame(sample_id = sample_ids, status = sample(c("A", "B"), length(sample_ids), replace = TRUE))

  task <- local_omnibus_task(list(notes = text_tbl), outcome)
  task <- grf_add_modality(task, "notes", "text")
  task <- grf_encode(task)

  encoded <- task$encodings$notes
  expect_s3_class(encoded, "tbl_df")
  expect_equal(nrow(encoded), length(sample_ids))
  expect_true(all(grepl("^text_", setdiff(names(encoded), "..rowid.."))))

  limited_task <- local_omnibus_task(list(notes = text_tbl), outcome)
  limited_task <- grf_add_modality(limited_task, "notes", grf_encoder_text(max_features = 3, term_count_min = 1))
  limited_task <- grf_encode(limited_task)
  limited <- limited_task$encodings$notes
  expect_lte(ncol(limited), 4)
})
