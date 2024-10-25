test_that("Initializing an empty summary_stats works", {
  expect_warning(as_summary_stats(data.table::data.table(), build = "b37")) |>
    expect_warning() |>
    expect_warning()
})


summary_stats <- data.table::data.table(
  variant_id = "chr1_2_A_T_b38",
  chr = "chr1",
  pos = 2,
  ref = "A",
  alt = "T",
  pval = .01,
  log10p = 2,
  effect = 1,
  effect_se = -1/qnorm(.01/2),
  zscore = -qnorm(.01/2),
  odds_ratio = 2
)
test_that("Standard initialization works with different column names", {
  expect_equal(as_summary_stats(summary_stats)[
    ,
    mget(colnames(summary_stats))
  ],
               summary_stats)
  column_names <- lapply(summary_stats_column_names,
                         sample,
                         size = 20,
                         replace = TRUE) |>
    data.table::as.data.table()

  apply(column_names, 1, \(names) {
    expect_equal(as_summary_stats(setNames(summary_stats, names))[
      ,
      mget(colnames(summary_stats))
    ],
                 summary_stats)
  })
})
