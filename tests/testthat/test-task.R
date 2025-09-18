test_that("grf_task aligns MultiAssayExperiment assays and outcome", {
  set.seed(123)
  ids <- paste0("S", 1:5)
  outcome <- data.frame(
    sample_id = sample(ids),
    status = factor(sample(c("A", "B"), length(ids), replace = TRUE))
  )

  mod1 <- matrix(rnorm(length(ids) * 3), nrow = 3, ncol = length(ids),
                 dimnames = list(paste0("g", 1:3), ids))
  mod2 <- matrix(rnorm(length(ids) * 2), nrow = 2, ncol = length(ids),
                 dimnames = list(paste0("p", 1:2), ids))

  se1 <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = mod1))
  se2 <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = mod2))
  mae <- MultiAssayExperiment::MultiAssayExperiment(
    experiments = MultiAssayExperiment::ExperimentList(list(omics = se1, proteomics = se2)),
    colData = S4Vectors::DataFrame(data.frame(dummy = rnorm(length(ids)), row.names = ids))
  )

  task <- grf_task(mae, outcome)

  expect_s3_class(task, "grf_task")
  aligned_ids <- rownames(MultiAssayExperiment::colData(task$mae))
  expect_equal(aligned_ids, outcome$sample_id)
  expect_equal(task$modalities, c("omics", "proteomics"))
  assays <- MultiAssayExperiment::experiments(task$mae)
  expect_equal(colnames(SummarizedExperiment::assay(assays[[1]])), aligned_ids)
  expect_equal(colnames(SummarizedExperiment::assay(assays[[2]])), aligned_ids)
  expect_equal(MultiAssayExperiment::colData(task$mae)[[task$id_col]], aligned_ids)
})

test_that("grf_task accepts lists of modalities", {
  ids <- paste0("P", 1:4)
  modalities <- list(
    omics = data.frame(sample_id = ids, f1 = rnorm(4), f2 = rnorm(4)),
    clinical = data.frame(sample_id = ids, age = 30:33, score = rnorm(4))
  )
  outcome <- data.frame(sample_id = ids, response = rnorm(4))

  task <- grf_task(modalities, outcome)

  expect_s3_class(task, "grf_task")
  assays <- MultiAssayExperiment::experiments(task$mae)
  expect_equal(colnames(SummarizedExperiment::assay(assays[["omics"]])), ids)
  expect_equal(colnames(SummarizedExperiment::assay(assays[["clinical"]])), ids)
  expect_equal(rownames(MultiAssayExperiment::colData(task$mae)), ids)
})

test_that("grf_task restricts to intersecting samples", {
  ids <- paste0("X", 1:4)
  modalities <- list(
    omics = data.frame(sample_id = ids, f1 = rnorm(4)),
    imaging = data.frame(sample_id = ids[1:3], pix = rnorm(3))
  )
  outcome <- data.frame(sample_id = ids, status = 1:4)

  task <- grf_task(modalities, outcome)

  expect_equal(rownames(MultiAssayExperiment::colData(task$mae)), ids[1:3])
  assays <- MultiAssayExperiment::experiments(task$mae)
  expect_equal(colnames(SummarizedExperiment::assay(assays[["omics"]])), ids[1:3])
  expect_equal(colnames(SummarizedExperiment::assay(assays[["imaging"]])), ids[1:3])
})

test_that("print.grf_task summarises the task", {
  ids <- paste0("M", 1:3)
  modalities <- list(
    omics = data.frame(sample_id = ids, f1 = 1:3, f2 = 4:6)
  )
  outcome <- data.frame(sample_id = ids, label = letters[1:3])

  task <- grf_task(modalities, outcome)

  expect_output(print(task), "Modalities \\(1\\): omics")
  expect_output(print(task), "Outcome columns: label")
})

test_that("built-in encoders produce tibbles with sample ids", {
  set.seed(99)
  ids <- sprintf("S%02d", 1:3)
  modalities <- list(
    omics = matrix(rnorm(length(ids) * 4), nrow = 4, dimnames = list(paste0("g", 1:4), ids)),
    clinical = data.frame(sample_id = ids, age = c(40, 55, 63), score = rnorm(length(ids))),
    notes = data.frame(sample_id = ids, text = c("tumor growth", "benign lesion", "tumor regression"))
  )
  outcome <- data.frame(sample_id = ids, status = factor(c("A", "B", "A")))

  task <- grf_task(modalities, outcome)
  task <- grf_add_modality(task, "omics", "numeric")
  task <- grf_add_modality(task, "clinical", "tabular")
  task <- grf_add_modality(task, "notes", "text")

  task <- grf_encode(task)

  expect_true(tibble::is_tibble(task$encodings$omics))
  expect_true(tibble::is_tibble(task$encodings$clinical))
  expect_true(tibble::is_tibble(task$encodings$notes))

  expect_setequal(task$encodings$omics$`..rowid..`, ids)
  expect_setequal(task$encodings$clinical$`..rowid..`, ids)
  expect_setequal(task$encodings$notes$`..rowid..`, ids)

  expect_gt(ncol(task$encodings$omics), 1)
  expect_true(ncol(task$encodings$clinical) >= 1)
  expect_true(ncol(task$encodings$notes) >= 1)
})
