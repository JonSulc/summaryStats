test_that("Initialization works", {
  expect_error(
    new_genomic_ranges()
  )
  expect_error(
    new_genomic_ranges(chromosome = "chr1")
  )
  expect_equal(
    new_genomic_ranges(
      chr = "chr1",
      start = 123,
      end = 234
    ),
    structure(
      data.table(
        chr = "chr1",
        start = 123,
        end = 234
      ),
      class = c("genomic_ranges", "data.table", "data.frame")
    )
  )
})

test_that("Validation works", {
  expect_no_error(
    validate_genomic_ranges(
      data.table::data.table(chr = NA_character_,
                             start = NA_integer_,
                             end = NA_integer_) |>
        as_genomic_ranges()
    )
  )
  expect_no_error(
    validate_genomic_ranges(
      data.table::data.table(chr = NA_character_,
                             pos = NA_integer_) |>
        as_genomic_ranges()
    )
  )
  expect_no_error(
    validate_genomic_ranges(
      data.table::data.table(chr = c("chr1", "chr2"),
                             start = c(123L, 1L),
                             end = c(234L, 2L)) |>
        as_genomic_ranges()
    )
  )
  expect_error(
    validate_genomic_ranges(
      data.table::data.table(chr = c("chr1", "chr2"),
                             start = c(123, 1),
                             end = c(234, 2)) |>
        as_genomic_ranges()
    ),
    "positive integers"
  )
  expect_error(
    validate_genomic_ranges(
      data.table::data.table(chr = c("chr1", "chr2"),
                             start = c(123L, 0L),
                             end = c(234L, 2L)) |>
        as_genomic_ranges()
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
        as_genomic_ranges()
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

test_that("Position is properly converted", {
  expect_error(assign_correct_start_end(1))
  expect_error(
    data.table::data.table(chr = "1") |>
      assign_correct_start_end()
  )
  expect_equal(
    data.table::data.table(start = "123", end = "234") |>
      assign_correct_start_end(),
    data.table::data.table(start = 123L, end = 234L)
  )
  expect_equal(
    data.table::data.table(start = 234, end = 123) |>
      assign_correct_start_end() |>
      suppressWarnings(),
    data.table::data.table(start = 123L, end = 234L)
  )
  expect_equal(
    data.table::data.table(start = "234", end = "123") |>
      assign_correct_start_end() |>
      suppressWarnings(),
    data.table::data.table(start = 123L, end = 234L)
  )
  expect_equal(
    data.table::data.table(
      start = c("234", 1),
      end   = c("123", 2)
    ) |>
      assign_correct_start_end() |>
      suppressWarnings(),
    data.table::data.table(start = c(123L, 1), end = c(234L, 2))
  )
})
