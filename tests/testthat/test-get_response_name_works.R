test_that("get_response_name returns the response", {
  expect_equal(get_response_name(a ~ b), "a")
})

test_that("get_response_name marks error without response in formula", {
  expect_error(get_response_name(~ b))
})