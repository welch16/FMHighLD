test_that("expect error when the input is not a formula", {

  expect_error(has_intercept("y ~ z"))

})

test_that("has_intercept works", {

  expect_identical(has_intercept(y ~ 1 + z), TRUE)
  expect_identical(has_intercept(y ~ z), TRUE)
  expect_identical(has_intercept(y ~ 0 + z), FALSE)

})
