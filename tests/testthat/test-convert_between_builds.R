library(quickcheck)

test_that("Build conversion works", {
  quickcheck::for_all(
    summary_stats = data_frame_(
      chr = integer_bounded(1, 22),
      pos = integer_bounded(1, 1e6),
      ref = character_letters(),
      alt = character_letters(),
      effect = quickcheck::double_(),
      pval = double_bounded(1e-10, 1)
    ),
    property = function(summary_stats) {
      previous <- as_summary_stats(summary_stats, build = "b37")
      capture.output(
        summary_stats <- as_summary_stats(summary_stats, build = "b37") |>
          convert_to_build("b38") |>
          convert_to_build("b37"),
        file = nullfile()
      )

      expect_equal(summary_stats[, .SD, .SDcols = c("variant_id", "chr", "pos", "ref", "alt")],
                   previous[, .SD, .SDcols = c("variant_id", "chr", "pos", "ref", "alt")])
    }
  )
})

test_that("Converting with NAs works", {
  expect_no_error(
    data.table::data.table(
      chr = NA,
      pos = NA,
      ref = "A",
      alt = "C",
      effect = 1,
      pval = 1
    ) |>
      as_summary_stats(build = "b38") |>
      convert_to_build("b37")
  )
  expect_no_error(
    data.table::data.table(
      chr = c(NA, "chr1"),
      pos = c(NA, 1),
      ref = "A",
      alt = "C",
      effect = 1,
      pval = 1
    ) |>
      as_summary_stats(build = "b38") |>
      convert_to_build("b37") |>
      suppressWarnings()
  )
})

test_that("Missing SNPs in build conversion are handled properly", {
  # This SNP does not have an equivalent in b37, should raise warning
  summary_stats <- data.table(
    variant_id = "chr5_49861645_A_G_b38",
    effect = 1,
    pval = .5
  ) |>
    as_summary_stats()
  expect_warning(capture.output(convert_to_build(summary_stats, "b37")),
                 "no mapped pos")
})
