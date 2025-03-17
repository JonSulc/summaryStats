test_that("Simulated stats initialize properly", {
  fdt <- fake_dt(
    nrows = 100,
    chr = 3,
    start = 42,
    end = 421
  )
  expect_equal(nrow(fdt), 100)
  expect_true(all(fdt[
    ,
    chr == "chr3" &
      42 <= pos & pos <= 421 &
      ref %chin% c("A", "C", "G", "T") & alt %chin% c("A", "C", "G", "T") &
      0 <= pval & pval <= 1 &
      is.numeric(effect)
  ]))
})
