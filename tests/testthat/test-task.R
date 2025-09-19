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

test_that("grf_fuse builds early, late, and hybrid structures", {
  ids <- sprintf("S%02d", 1:3)
  modalities <- list(
    omics = matrix(rnorm(length(ids) * 5), nrow = 5, dimnames = list(paste0("g", 1:5), ids)),
    clinical = data.frame(sample_id = ids, age = c(41, 59, 66), score = rnorm(length(ids))),
    notes = data.frame(sample_id = ids, text = c("tumor growth rapid", "benign lesion noted", "tumor shrinking"))
  )
  outcome <- data.frame(sample_id = ids, status = factor(c("A", "B", "A")))

  task <- grf_task(modalities, outcome)
  task <- grf_add_modality(task, "omics", "numeric")
  task <- grf_add_modality(task, "clinical", "tabular")
  task <- grf_add_modality(task, "notes", "text")
  task <- grf_encode(task)

  early <- grf_fuse(task, "early")
  expect_equal(early$fusion$mode, "early")
  expect_true(tibble::is_tibble(early$fusion$X))
  expect_equal(nrow(early$fusion$X), length(ids))
  expect_setequal(early$fusion$X$`..rowid..`, ids)
  expect_true(all(grepl("^(omics|clinical|notes)__", setdiff(names(early$fusion$X), "..rowid.."))))
  expect_null(early$fusion$per_modality)

  late <- grf_fuse(task, "late")
  expect_equal(late$fusion$mode, "late")
  expect_null(late$fusion$X)
  expect_length(late$fusion$per_modality, 3)
  expect_named(late$fusion$per_modality, c("omics", "clinical", "notes"))
  purrr::iwalk(late$fusion$per_modality, function(tbl, nm){
    expect_true(tibble::is_tibble(tbl))
    expect_equal(nrow(tbl), length(ids))
    expect_setequal(tbl$`..rowid..`, ids)
    expect_true("..rowid.." %in% names(tbl))
  })

  hybrid <- grf_fuse(task, "hybrid")
  expect_equal(hybrid$fusion$mode, "hybrid")
  expect_true(tibble::is_tibble(hybrid$fusion$X))
  expect_equal(hybrid$fusion$X, early$fusion$X)
  expect_equal(hybrid$fusion$per_modality, late$fusion$per_modality)
})

test_that("grf_fuse retains intersection of encoded sample ids", {
  ids <- sprintf("S%02d", 1:4)
  modalities <- list(
    omics = matrix(rnorm(length(ids) * 3), nrow = 3, dimnames = list(paste0("g", 1:3), ids)),
    clinical = data.frame(sample_id = ids, age = seq_len(length(ids)))
  )
  outcome <- data.frame(sample_id = ids, status = rnorm(length(ids)))

  drop_first_encoder <- function(se, id_col) {
    sample_ids <- colnames(SummarizedExperiment::assay(se))
    keep_ids <- sample_ids[1:3]
    tibble::tibble(`..rowid..` = keep_ids, feature = seq_along(keep_ids))
  }

  drop_last_encoder <- function(se, id_col) {
    sample_ids <- colnames(SummarizedExperiment::assay(se))
    keep_ids <- sample_ids[2:4]
    tibble::tibble(`..rowid..` = keep_ids, value = seq_along(keep_ids))
  }

  task <- grf_task(modalities, outcome)
  task <- grf_add_modality(task, "omics", drop_first_encoder)
  task <- grf_add_modality(task, "clinical", drop_last_encoder)
  task <- grf_encode(task)

  fused <- grf_fuse(task, "hybrid")
  expect_equal(fused$fusion$X$`..rowid..`, ids[2:3])
  expect_equal(fused$fusion$per_modality$omics$`..rowid..`, ids[2:3])
  expect_equal(fused$fusion$per_modality$clinical$`..rowid..`, ids[2:3])
  expect_equal(nrow(fused$fusion$X), 2)
})
