test_that("are_only_positive_integers throws errors when necessary", {
  expect_error(
    are_only_positive_integers(1, "test")
  )
  expect_error(
    are_only_positive_integers(data.frame(test = 1), "test")
  )
  expect_error(
    are_only_positive_integers(data.table::data.table(wrong = 1.), "wrong"),
    "integers"
  )
  expect_error(
    are_only_positive_integers(data.table::data.table(wrong = -1L), "wrong"),
    "positive"
  )
  expect_error(
    are_only_positive_integers(data.table::data.table(wrong = 0L), "wrong"),
    "positive"
  )
  expect_error(
    data.table::data.table(wrong = 0L, correct = 1L) |>
      are_only_positive_integers(c("wrong", "correct")),
    "^(?!.*correct).*wrong.*$",
    perl = TRUE
  )
  expect_no_error(
    are_only_positive_integers(data.table::data.table(correct = 1L), "correct")
  )
  expect_equal(
    are_only_positive_integers(data.table::data.table(correct = 1L), "correct"),
    data.table::data.table(correct = 1L)
  )
})

test_that("make_positive_integer works", {
  expect_equal(
    make_positive_integer(1),
    1L
  )
  expect_equal(
    make_positive_integer(-1),
    NA_integer_
  )
  expect_equal(
    make_positive_integer("1"),
    1L
  )
  expect_equal(
    make_positive_integer("test") |>
      suppressWarnings(),
    NA_integer_
  )
  expect_equal(
    make_positive_integer(1:3),
    1L:3L
  )
  expect_equal(
    make_positive_integer(-1:1),
    c(NA_integer_, NA_integer_, 1L)
  )
  expect_equal(
    make_positive_integer(c("1", NA, 42)),
    c(1L, NA_integer_, 42L)
  )
})
