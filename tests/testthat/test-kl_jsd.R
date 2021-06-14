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

test_that("jsd gets zero for same vector", {

  p1 <- c(3, 2, 1)
  expect_equal(jsd(p1, p1), 0)

})

test_that("kl works when an entry is zero", {

  p1 <- 1:10
  p1[5] <- 0
  p2 <- rep(1, length(p1))

  # this test was motivated by the issue that is.na(sum(p1 * log(p1))) is TRUE

  expect_true(! is.na(kl(p1, p2)))

})

test_that("kl error when length <= 1", {

  expect_error(kl(0, 1))

})

test_that("kl error when the vectors are of diff length", {

  expect_error(kl(rnorm(3), rnorm(2)))

})

test_that("kl works when empty", {

  expect_equal(kl(rep(0, 2), rep(1, 2)), 0)

})

test_that("jsd works for two vectors", {

  p1 <- 1:5
  p1 <- p1 / sum(p1)
  q1 <- c(4, 3, 2, 1, 4)
  q1 <- q1 / sum(q1)

  m1 <- colMeans(rbind(p1, q1))

  expect_equivalent(jsd(p1, q1), mean(c(kl(p1, m1), kl(q1, m1))))

})
