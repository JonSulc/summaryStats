library(quickcheck)

test_that("Validation detects any invalid entries", {
  expect_error(
    validate_summary_stats(1),
    "inherits"
  )
  expect_error(
    validate_summary_stats(structure(1, class = "summary_stats"))
  )

  summary_stats <- fake_stats()
  expect_no_error(validate_summary_stats(summary_stats))
  expect_equal(validate_summary_stats(summary_stats), summary_stats)
  for (column in c("variant_id", "chr", "pos", "ref", "alt", "effect", "pval", "effect_se")) {
    expect_error(
      validate_summary_stats(summary_stats[, .SD, .SDcols = -column]),
      paste("missing:", column)
    )
  }

  expect_error(
    summary_stats[
      ,
      c(.SD,
        variant_id = "duplicates"),
      .SDcols = -"variant_id"
    ] |>
      validate_summary_stats(),
    "not unique"
  )

  expect_error(
    summary_stats[
      ,
      c(.SD,
        .(chr = 1)),
      .SDcols = -"chr"
    ] |>
      validate_summary_stats(),
    "chr"
  )

  expect_no_error(
    summary_stats[
      ,
      c(.SD,
        .(pos = 1.)),
      .SDcols = -"pos"
    ] |>
      validate_summary_stats()
  )
  expect_error(
    summary_stats[
      ,
      c(.SD,
        .(pos = -1)),
      .SDcols = -"pos"
    ] |>
      validate_summary_stats(),
    "Non-positive"
  )
  expect_error(
    summary_stats[
      ,
      c(.SD,
        .(pos = "here")),
      .SDcols = -"pos"
    ] |>
      validate_summary_stats(),
    "numeric"
  )
  expect_error(
    summary_stats[
      ,
      c(.SD,
        .(pos = pi)),
      .SDcols = -"pos"
    ] |>
      validate_summary_stats(),
    "Non-integer"
  )

  expect_error(
    summary_stats[
      ,
      c(.SD,
        .(ref = 1)),
      .SDcols = -"ref"
    ] |>
      validate_summary_stats(),
    "Invalid"
  )
  expect_error(
    summary_stats[
      ,
      c(.SD,
        .(ref = "Something")),
      .SDcols = -"ref"
    ] |>
      validate_summary_stats(),
    "Invalid"
  )
  expect_no_error(
    summary_stats[
      ,
      c(.SD,
        .(ref = "ACTGGTCA")),
      .SDcols = -"ref"
    ] |>
      validate_summary_stats()
  )

  expect_error(
    summary_stats[
      ,
      c(.SD,
        .(alt = 1)),
      .SDcols = -"alt"
    ] |>
      validate_summary_stats(),
    "Invalid"
  )
  expect_error(
    summary_stats[
      ,
      c(.SD,
        .(alt = "Something")),
      .SDcols = -"alt"
    ] |>
      validate_summary_stats(),
    "Invalid"
  )
  expect_no_error(
    summary_stats[
      ,
      c(.SD,
        .(alt = "ACTGGTCA")),
      .SDcols = -"alt"
    ] |>
      validate_summary_stats()
  )

  expect_error(
    summary_stats[
      ,
      c(.SD,
        .(alt = ref)),
      .SDcols = -"alt"
    ] |>
      validate_summary_stats(),
    "identical"
  )

  expect_error(
    summary_stats[
      ,
      c(.SD,
        .(effect = "NA")),
      .SDcols = -"effect"
    ] |>
      validate_summary_stats(),
    "Non-numeric"
  )

  expect_error(
    summary_stats[
      ,
      c(.SD,
        .(effect_se = "NA")),
      .SDcols = -"effect_se"
    ] |>
      validate_summary_stats(),
    "Non-numeric"
  )
  expect_error(
    summary_stats[
      ,
      c(.SD,
        .(effect_se = 0)),
      .SDcols = -"effect_se"
    ] |>
      validate_summary_stats(),
    "0"
  )
})

test_that("Regular initialization of summary_stats works", {
  dt <- data.table::data.table(
    chr = 1,
    pos = 123,
    ref = "A",
    alt = "C",
    effect = 1,
    pval = .05
  )
  summary_stats <- data.table::data.table(
    variant_id = "chr1_123_A_C_b38",
    chr = "chr1",
    pos = 123,
    ref = "A",
    alt = "C",
    effect = 1,
    pval = .05,
    effect_se = 1/abs(qnorm(.025)),
    chr_b38 = "chr1",
    pos_b38 = 123
  ) |>
    data.table::setattr("class", c("summary_stats", "data.table", "data.frame"))

  expect_equal(
    new_summary_stats(dt, build = "b38"),
    summary_stats
  )

  expect_equal(
    new_summary_stats(dt[, lapply(.SD, \(x) {
      if (is.character(x)) tolower(x) else x
    })], build = "b38"),
    summary_stats
  )

  summary_stats <- fake_stats()
  expect_equal(
    summary_stats[
      ,
      2*pnorm(-abs(effect/effect_se))
    ],
    summary_stats$pval
  )
  expect_equal(
    summary_stats[
      ,
      paste(chr, pos, ref, alt, "b38", sep = "_")
    ],
    summary_stats$variant_id
  )
})

test_that("Initializing empty summary_stats works", {
  expect_warning(
    new_summary_stats(data.table::data.table(),
                         build = "b37")
  )
})

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

test_that("Standard initialization works with different column names", {
  summary_stats <- fake_stats() |>
    new_summary_stats()
  for (i in 1:10) {
    new_names <- sapply(
      colnames(summary_stats),
      \(cname) sample(summary_stats_column_names[[cname]], 1)
    )
    expect_equal(
      setNames(summary_stats,
               new_names) |>
        new_summary_stats(),
      summary_stats
    )
  }
})
