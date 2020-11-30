test_that("jsd gets zero for same vector", {

  p1 <- c(3, 2, 1)
  expect_equal(jsd(p1, p1), 0)

})
