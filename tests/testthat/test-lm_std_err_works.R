test_that("lm std. error works", {

  n <- 50
  mydata <- tibble::tibble(x = rnorm(n))
  mydata$y <- 5 + 3 * mydata$x + rnorm(n, 0, 3)
  mod <- lm(y ~ x, data = mydata)

  expect_equal(lm_std_err(mod),
    coef(summary(mod))[, "Std. Error"])
  expect_vector(lm_std_err(mod))
})


test_that("error when lm is not linear model", {

  expect_error(lm_std_err(10))
  expect_error(lm_std_err(c("a", "b")))
})
