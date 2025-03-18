test_that("Getting column names works", {
  expect_equal(
    get_standard_colnames(letters[1:3],
                          setNames(letters[1:3],
                                   as.list(paste0("col", 1:3)))),
    setNames(as.list(letters[1:3]),
             paste0("col", 1:3))
  )

  for (i in 1:10) {
    bad_names <- sapply(
      sample(names(summary_stats_column_names),
             sample(seq_along(summary_stats_column_names), 1)),
      \(cname) sample(summary_stats_column_names[[cname]], 1)
    )
    standardized_colnames <- get_standard_colnames(bad_names)
    expect_true(
      all(names(standardized_colnames) == names(summary_stats_column_names))
    )
    expect_identical(
      setdiff(unlist(standardized_colnames), bad_names),
      character(0)
    )
    expect_true(
      all(
        sapply(names(unlist(standardized_colnames)),
               \(cname) standardized_colnames[[cname]] %in% summary_stats_column_names[[cname]])
      )
    )
  }
})

test_that("Parsing variant_id works", {
  expect_false(is_correct_variant_id_format(NULL))
  expect_true(is_correct_variant_id_format("chr1_123_A_G_b38"))
  expect_false(is_correct_variant_id_format("chr1_123_a_g_b38"))
  expect_false(is_correct_variant_id_format("chr1_123_N_T_b38"))
  expect_false(is_correct_variant_id_format("rs1234"))
})

test_that("Getting build colnames works", {
  expect_equal(
    get_build_colnames("b38", columns = "test"),
    "test_b38"
  )
  expect_equal(
    get_build_colnames("b38", columns = c("variant_id", "chr", "pos")),
    c("variant_id_b38", "chr_b38", "pos_b38")
  )
})

test_that("Checking build version works", {
  expect_true(has_build_version(fake_stats(), "b38"))
  expect_false(has_build_version(fake_stats(), "b37"))
  expect_true(has_build_version(fake_stats(build = "b37"), "b37"))
  expect_false(has_build_version(fake_stats(build = "b37"), "b38"))
})

test_that("Genome build can be inferred from summary statistics", {
  expect_equal(
    structure("test", build = "b38") |>
      get_build(),
    "b38"
  )
  expect_equal(
    get_build(fake_stats()),
    "b38"
  )
  expect_equal(
    fake_stats() |>
      data.table::setattr("build", NULL) |>
      get_build(),
    "b38"
  )
  expect_equal(
    fake_stats(build = "b37") |>
      data.table::setattr("build", NULL) |>
      get_build(),
    "b37"
  )
  expect_equal(
    fake_stats(0) |>
      data.table::setattr("build", NULL) |>
      get_build() |>
      suppressWarnings(),
    "b38"
  )
  expect_equal(
    fake_stats(0, build = "b37") |>
      data.table::setattr("build", NULL) |>
      get_build() |>
      suppressWarnings(),
    "b37"
  )
  expect_equal(
    fake_stats()[, -c("variant_id_b38", "chr_b38", "pos_b38")] |>
      data.table::setattr("build", NULL) |>
      get_build(),
    "b38"
  )
})

test_that("Position can be inferred from non-variant_id columns", {
  expect_error(
    data.table::data.table(chr = 1, ref = "a", alt = "c") |>
      infer_position(),
    "Unable to infer"
  )
  expect_equal(
    data.table::data.table(MarkerName = "chr1:123:a:c") |>
      infer_position(),
    123
  )
})

test_that("Parsing variant_id produces chr, pos, ref, alt", {
  expect_equal(
    parse_variant_id("chr1_123_A_C_b38"),
    data.table::data.table(chr = "chr1", pos = 123, ref = "A", alt = "C")
  )
  expect_error(
    parse_variant_id("chr1:123:A:C"),
    "Unable to parse"
  )
})

test_that("Generating variant_id from chr, pos, ref, alt works", {
  expect_equal(
    get_variant_id_from_chr_pos_ref_alt(
      data.table::data.table(
        chr = "chr1",
        pos = 123,
        ref = "A",
        alt = "C"
      ),
      build = "b38"
    ),
    "chr1_123_A_C_b38"
  )
  expect_equal(
    get_variant_id_from_chr_pos_ref_alt(
      data.table::data.table(
        chr = "chr1",
        pos = c(123, 234),
        ref = "A",
        alt = "C"
      ),
      build = "b38"
    ),
    c("chr1_123_A_C_b38", "chr1_234_A_C_b38")
  )
  expect_error(
    get_variant_id_from_chr_pos_ref_alt(
      data.table::data.table(
        chr = "chr1",
        pos = c(123, 234),
        ref = "A"
      ),
      build = "b38"
    ),
    "Missing.*: alt$"
  )
})
