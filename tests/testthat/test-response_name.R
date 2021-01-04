test_that("get right response name", {
  expect_equal(get_response_name(y ~ x), "y")
})

test_that("get error without response name", {
  expect_error(get_response_name(~ x))
})
