library(quickcheck)

test_that("Initializing an empty summary_stats works", {
  expect_warning(as_summary_stats(data.table::data.table(), build = "b37")) |>
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


test_that("Getting column names works", {
  all_colnames <- setNames(letters, paste0("V", 1:26)) |>
    as.list()

  for_all(
    indices = integer_bounded(1, 26, len = 1:10),
    property = function(indices) {
      summary_stats <- matrix(0, ncol = length(indices), nrow = 5) |>
        data.table::as.data.table() |>
        setNames(letters[indices])

      columns <- get_colnames(summary_stats, all_colnames)

      expect_equal(columns[indices],
                   as.list(setNames(letters[indices], paste0("V", indices))))
      expect_true(all(sapply(columns[-indices], length) == 0))
    }
  )
})
