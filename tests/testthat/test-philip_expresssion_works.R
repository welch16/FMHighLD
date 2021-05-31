
test_that("philips works", {

  vec <- 1:5

  expect_equal(philips(rep(0, 5)), 0)
  expect_equal(philips(vec), sqrt(sum(vec^2) / max(1 + vec)))

})

test_that("error when non numeric", {
  expect_error(philips(c("a", "b")))
})

test_that("work with logical vectors", {

  expect_equal(philips(c(1, 0, 1, 1)),
    philips(c(TRUE, FALSE, TRUE, TRUE)))

})