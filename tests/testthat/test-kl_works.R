test_that("kl gets zero for same vector", {

  p1 <- c(3, 2, 1)

  expect_equal(kl(p1, p1), 0)
})
