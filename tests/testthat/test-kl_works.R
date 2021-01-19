test_that("kl gets zero for same vector", {

  p1 <- c(3, 2, 1)

  expect_equal(kl(p1, p1), 0)
})

test_that("kl works for two vectors", {

  p1 <- 1:5
  p1 <- p1 / sum(p1)
  q1 <- c(4, 3, 2, 1, 4)
  q1 <- q1 / sum(q1)

  expect_equivalent(kl(p1, q1), sum(p1 * log(p1 / q1)))

})
