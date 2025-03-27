test_that("'chr' validation works", {
  expect_error(validate_chr.genomic_positions(2))
  expect_error(validate_chr.genomic_positions(
    data.table(chr = 1)
  ),
  "character")
  expect_error(validate_chr.genomic_positions(
    data.table(chr = "1")
  ),
  "start with 'chr'")
  expect_no_error(validate_chr.genomic_positions(
    data.table(chr = "chr1")
  ))
})

test_that("is_any_non_positive works", {
  expect_true(
    is_any_non_positive(-1)
  )
  expect_false(
    is_any_non_positive(1)
  )
  expect_true(
    is_any_non_positive(c(-1, NA))
  )
  expect_false(
    is_any_non_positive(c(1, NA))
  )
  expect_false(
    is_any_non_positive(NA)
  )
})

test_that("'pos' validation works", {
  expect_error(
    validate_pos.genomic_positions(1),
    "NULL"
  )
  expect_error(
    validate_pos.genomic_positions(
      structure(1, pos = "pos")
    ),
    "Columns"
  )
  expect_error(
    validate_pos.genomic_positions(
      data.frame(pos = "not numeric") |>
        structure(pos = "pos"),
      "integer"
    )
  )
  expect_error(
    validate_pos.genomic_positions(
      data.frame(pos = 0) |>
        structure(pos = "pos"),
      "positive"
    )
  )
  expect_no_error(
    validate_pos.genomic_positions(
      data.table::data.table(pos = 1L) |>
        structure(pos = "pos")
    )
  )
  expect_equal(
    validate_pos.genomic_positions(
      data.table::data.table(pos = 1L) |>
        structure(pos = "pos")
    ),
    data.table::data.table(pos = 1L) |>
      structure(pos = "pos")
  )
})

test_that("Initialization works", {
  expect_error(
    new_genomic_positions(),
    "missing"
  )
  expect_error(
    new_genomic_positions(data = "chr1"),
    "data.frame"
  )
  expect_equal(
    new_genomic_positions(
      data.frame(chr   = "chr1",
                 start = 123,
                 end   = 234)
    ),
    structure(
      data.table::data.table(
        chr   = "chr1",
        start = 123L,
        end   = 234L
      ),
      pos = c("start", "end"),
      build = "b38",
      class = c("genomic_positions", "data.table", "data.frame")
    )
  )

  expect_warning(
    new_genomic_positions(
      data.frame(chr   = "chr1",
                 start = -123,
                 end   = 234)
    ),
    "not positive"
  )
  expect_equal(
    new_genomic_positions(
      data.frame(chr   = "chr1",
                 start = -123,
                 end   = 234)
    ) |>
      suppressWarnings(),
    structure(
      data.table::data.table(
        chr   = "chr1",
        start = NA_integer_,
        end   = 234L
      ),
      pos = c("start", "end"),
      build = "b38",
      class = c("genomic_positions", "data.table", "data.frame")
    )
  )
})
