test_that("select kth works when length is 1", {

  res <- 1

  expect_equal(select_kth_random(res, .95, 2), 1)
})

test_that("select kth works when length > 1", {

  set.seed(12312)
  res <- sort(abs(rnorm(10)))
  expect_true(select_kth_random(res, .95, 2) %in% 1:2)

})
