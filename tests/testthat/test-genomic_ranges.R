test_that("Initialization works", {
  expect_error(
    new_genomic_ranges()
  )
  expect_error(
    new_genomic_ranges(chromosome = "chr1")
  )
  expect_equal(
    new_genomic_ranges(
      data.table::data.table(
        chr   = "chr1",
        start = 123,
        end   = 234
      )
    ),
    structure(
      data.table::data.table(
        chr   = "chr1",
        start = 123,
        end   = 234
      ),
      pos   = c("start", "end"),
      build = "b38",
      class = c("genomic_ranges", "genomic_positions",
                "data.table", "data.frame")
    )
  )
})

test_that("Validation works", {
  expect_no_error(
    validate_genomic_ranges(
      data.table::data.table(chr = NA_character_,
                             start = NA_integer_,
                             end = NA_integer_) |>
        structure(
          class = c("genomic_ranges", "genomic_positions",
                    "data.table", "data.frame"),
          build = "b38",
          pos   = c("start", "end")
        )
    )
  )
  expect_error(
    validate_genomic_ranges(
      data.table::data.table(chr = NA_character_,
                             pos = NA_integer_) |>
        structure(
          class = c("genomic_ranges", "genomic_positions",
                    "data.table", "data.frame"),
          build = "b38",
          pos   = c("start", "end")
        )
    ),
    "start.*end"
  )
  expect_no_error(
    validate_genomic_ranges(
      data.table::data.table(chr   = c("chr1", "chr2"),
                             start = c(123L, 1L),
                             end   = c(234L, 2L)) |>
        structure(
          class = c("genomic_ranges", "genomic_positions",
                    "data.table", "data.frame"),
          build = "b38",
          pos   = c("start", "end")
        )
    )
  )
  expect_error(
    validate_genomic_ranges(
      data.table::data.table(chr   = c("chr1", "chr2"),
                             start = c(123, 1),
                             end   = c(234, 2)) |>
        structure(
          class = c("genomic_ranges", "genomic_positions",
                    "data.table", "data.frame"),
          build = "b38",
          pos   = c("start", "end")
        )
    ),
    "positive integers"
  )
  expect_error(
    validate_genomic_ranges(
      data.table::data.table(chr = c("chr1", "chr2"),
                             start = c(123L, 0L),
                             end = c(234L, 2L)) |>
        structure(
          class = c("genomic_ranges", "genomic_positions",
                    "data.table", "data.frame"),
          build = "b38",
          pos   = c("start", "end")
        )
    ),
    "positive integers"
  )
  expect_error(
    validate_genomic_ranges(1),
    "columns"
  )
  expect_error(
    validate_genomic_ranges(
      data.table::data.table(chr = c("chr1", "chr2"))
    ),
    "columns"
  )
  expect_error(
    validate_genomic_ranges(
      data.table::data.table(end = c(234, 2))
    ),
    "columns"
  )
  expect_error(
    validate_genomic_ranges(
      data.table::data.table(chr = 1:2,
                             start = c(123, 1),
                             end = c(234, 2)) |>
        as_genomic_ranges()
    ),
    "character"
  )
  expect_error(
    validate_genomic_ranges(
      data.table::data.table(chr = c("1", "chr2"),
                             start = c(123, 1),
                             end = c(234, 2)) |>
        as_genomic_ranges()
    ),
    "'chr'"
  )
  expect_error(
    validate_genomic_ranges(
      data.table::data.table(chr = c("chr1", "chr2"),
                             start = c(1234L, 1L),
                             end = c(234L, 2L)) |>
        structure(
          class = c("genomic_ranges", "genomic_positions",
                    "data.table", "data.frame"),
          build = "b38",
          pos   = c("start", "end")
        )
    ),
    "lower"
  )
})

test_that("Assigning chr works with genomic ranges", {
  expect_equal(
    data.table::data.table(
      chr = "1"
    ) |>
      assign_correct_chr(),
    data.table::data.table(
      chr = "chr1"
    )
  )
  expect_equal(
    data.table::data.table(
      chr = 1
    ) |>
      assign_correct_chr(),
    data.table::data.table(
      chr = "chr1"
    )
  )
})
