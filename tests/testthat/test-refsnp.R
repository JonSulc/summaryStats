test_that("Reading JSON file works", {
  expect_no_error(
    read_json_file(
      "~/rcp_storage/common/Users/sulc/data/dbsnp/refsnp-chr22.json.bz2",
      nlines = 10
    )
  )
  expect_equal(
    read_json_file(
      "~/rcp_storage/common/Users/sulc/data/dbsnp/refsnp-chr22.json.bz2",
      nlines = 10
    )[0],
    data.table::data.table(
      build = character(0),
      chr   = character(0),
      pos   = numeric(0),
      ref   = character(0),
      alt   = character(0)
    )
  )
})

refsnp_test <- data.table::data.table(
  build = rep(c("b37", "b38"), 10),
  chr = "chr22",
  pos = 1:20,
  ref = "A",
  alt = "C"
)
genomic_ranges <- data.table::data.table(
  chr = "chr22",
  pos = 2*5:15,
  ref = "A",
  alt = "C"
)
test_that("Inferring build version works", {
  expect_equal(
    infer_build(genomic_ranges, refsnp_test) |>
      withr::with_output_sink(new = "/dev/null"),
    "b38"
  )
  expect_equal(
    infer_build(genomic_ranges[1:5], refsnp_test) |>
      withr::with_output_sink(new = "/dev/null"),
    "b38"
  )
  expect_equal(
    infer_build(genomic_ranges[6:11], refsnp_test) |>
      withr::with_output_sink(new = "/dev/null"),
    "b36"
  )
  expect_equal(
    infer_build(genomic_ranges, refsnp_test, TRUE) |>
      withr::with_output_sink(new = "/dev/null"),
    data.table::data.table(
      build = c("b38", "b36"),
      build_match = c(6:5/11)
    )
  )
  expect_equal(
    infer_build(genomic_ranges[1:5], refsnp_test, TRUE) |>
      withr::with_output_sink(new = "/dev/null"),
    data.table::data.table(
      build = c("b38"),
      build_match = 1
    )
  )
  expect_equal(
    infer_build(genomic_ranges[6:11], refsnp_test, TRUE) |>
      withr::with_output_sink(new = "/dev/null"),
    data.table::data.table(
      build = c("b38", "b36"),
      build_match = c(1, 5)/6
    )
  )
})

test_that("Builds inferred can sum up to more than 100%", {
  expect_equal(
    infer_build(genomic_ranges[1:5],
                rbind(refsnp_test,
                      refsnp_test[
                        ,
                        c(.(build = data.table::fifelse(build == "b38", "b37", "b38")),
                          .SD),
                        .SDcols = -"build"
                      ]),
                TRUE) |>
      withr::with_output_sink(new = "/dev/null"),
    data.table::data.table(
      build = c("b38", "b37"),
      build_match = 1
    )
  )
})
