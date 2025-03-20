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

  expect_error(
    summary_stats[
      ,
      c(.SD,
        .(pos = 1.)),
      .SDcols = -"pos"
    ] |>
      validate_summary_stats(),
    "integer"
  )
  expect_error(
    summary_stats[
      ,
      c(.SD,
        .(pos = "here")),
      .SDcols = -"pos"
    ] |>
      validate_summary_stats(),
    "integer"
  )
  expect_error(
    summary_stats[
      ,
      c(.SD,
        .(pos = pi)),
      .SDcols = -"pos"
    ] |>
      validate_summary_stats(),
    "integer"
  )
  expect_error(
    summary_stats[
      ,
      c(.SD,
        .(pos = -1L)),
      .SDcols = -"pos"
    ] |>
      validate_summary_stats(),
    "positive"
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
    variant_id_b38 = "chr1_123_A_C_b38",
    chr_b38 = "chr1",
    pos_b38 = 123
  ) |>
    data.table::setattr("class", c("summary_stats", "data.table", "data.frame")) |>
    data.table::setattr("build", "b38")

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
  expect_equal(
    empty_summary_stats(),
    set_as_summary_stats(
      data.table::data.table(
        variant_id = character(0),
        chr = character(0),
        pos = integer(0),
        ref = character(0),
        alt = character(0),
        effect = numeric(0),
        effect_se = numeric(0),
        pval = numeric(0)
      )
    )
  )
  expect_equal(
    empty_summary_stats(build = "b38"),
    set_as_summary_stats(
      data.table::data.table(
        variant_id = character(0),
        chr = character(0),
        pos = integer(0),
        ref = character(0),
        alt = character(0),
        effect = numeric(0),
        effect_se = numeric(0),
        pval = numeric(0),
        variant_id_b38 = character(0),
        chr_b38 = character(0),
        pos_b38 = integer(0)
      )
    ) |>
      data.table::setattr("build", "b38")
  )
  expect_warning(
    new_summary_stats(data.table::data.table(),
                         build = "b37")
  )
})

test_that("Standard initialization works with different column names", {
  summary_stats <- fake_stats() |>
    new_summary_stats()
  summary_stats <- summary_stats[, .SD, .SDcol = !names(summary_stats) %like% "_b38$"]
  for (i in 1:10) {
    new_names <- sapply(
      colnames(summary_stats),
      \(cname) sample(summary_stats_column_names[[cname]], 1)
    )
    expect_equal(
      setNames(summary_stats,
               new_names) |>
        new_summary_stats(),
      summary_stats |>
        new_summary_stats()
    )
  }
})

test_that("Inferring missing ref/alt alleles works", {
  expect_equal(
    data.table::data.table(
      ref = "A",
      alt = "C"
    ) |>
      infer_ref_alt_alleles(),
    data.table::data.table(
      ref = "A",
      alt = "C"
    )
  )
  expect_equal(
    data.table::data.table(
      a1 = c("A", "A"),
      a2 = c("C", "C"),
      alt = c("A", "C")
    ) |>
      infer_ref_alt_alleles(),
    data.table::data.table(
      ref = c("C", "A"),
      alt = c("A", "C")
    )
  )
  expect_equal(
    data.table::data.table(
      a1 = c("A", "A"),
      a2 = c("C", "C"),
      alt = c("C", "C")
    ) |>
      infer_ref_alt_alleles(),
    data.table::data.table(
      ref = c("A", "A"),
      alt = c("C", "C")
    )
  )

  expect_equal(
    data.table::data.table(
      a1 = "A",
      a2 = "C"
    ) |>
      infer_ref_alt_alleles() |>
      suppressWarnings(),
    data.table::data.table(
      ref = "A",
      alt = "C",
      assumed_ea_col = "a2"
    )
  )

  expect_warning(
    data.table::data.table(
      a1 = "A",
      a2 = "C"
    ) |>
      infer_ref_alt_alleles(),
    "assuming"
  )
  expect_error(
    data.table::data.table(
      al1 = "A",
      al2 = "C"
    ) |>
      infer_ref_alt_alleles(),
    "Failed"
  )
})

test_that("Completing summary stats with inferred alleles works", {
  expect_equal(
    data.table::data.table(
      a1 = "A",
      a2 = "C"
    ) |>
      infer_and_assign_ref_alt_alleles() |>
      suppressWarnings(),
    data.table::data.table(
      a1 = "A",
      a2 = "C",
      ref = "A",
      alt = "C",
      assumed_ea_col = "a2"
    )
  )
})

test_that("Any missing information is inferred", {
  summary_stats <- fake_stats()
  for (cname in c("variant_id", "chr", "pos", "ref", "alt")) {
    expect_equal(
      summary_stats[, .SD, .SDcols = -cname] |>
        new_summary_stats(build = "b38"),
      summary_stats
    )
  }
  expect_equal(
    summary_stats[, -c("chr", "pos", "ref", "alt")] |>
      new_summary_stats(build = "b38"),
    summary_stats
  )
})
