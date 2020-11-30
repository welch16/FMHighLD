test_that("kl works", {

  p1 <- c(3, 2, 1)

  expect_equal(kl(p1, p1), 0)
})
