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
