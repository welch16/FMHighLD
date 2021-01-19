test_that("jsd gets zero for same vector", {

  p1 <- c(3, 2, 1)
  expect_equal(jsd(p1, p1), 0)

})

test_that("jsd works for two vectors", {

  p1 <- 1:5
  p1 <- p1 / sum(p1)
  q1 <- c(4, 3, 2, 1, 4)
  q1 <- q1 / sum(q1)

  m1 <- colMeans(rbind(p1, q1))

  expect_equivalent(jsd(p1, q1), mean(c(kl(p1, m1), kl(q1, m1))))

})
