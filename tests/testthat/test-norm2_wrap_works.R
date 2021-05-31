test_that("norm2_wrap works", {

  vec <- as.matrix(1:10)

  expect_equal(norm2_wrap(rep(0, 10)), 0)
  expect_equal(norm2_wrap(rep(1, 5)), sqrt(5))
  expect_equal(norm2_wrap(vec), sqrt(sum(vec^2)))

})

test_that("error non-numeric vector", {


  expect_error(norm2_wrap(c("a", "b", "c")))


})