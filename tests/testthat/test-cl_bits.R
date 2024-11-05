test_that("Base file reading works", {
  expect_true(
    grepl("^zcat ", cl_read_dbgap_file_template("test", "mvp"))
  )
  expect_true(
    grepl("^awk ", cl_read_dbgap_file_template("test", "charge"))
  )
})
