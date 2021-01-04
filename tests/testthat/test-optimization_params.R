
fm1 <- FMParam()
fm2 <- FMParam(error_bound = 1e-3, max_iter = 200, min_tol = 1e-2)

test_that("Get correct error bound", {
  expect_equivalent(error_bound(fm1), 1e-4)
  expect_equivalent(error_bound(fm2), 1e-3)
})

test_that("Get correct max iter", {
  expect_equivalent(max_iter(fm1), 100)
  expect_equivalent(max_iter(fm2), 200)
})

test_that("Get correct min tol", {
  expect_equivalent(min_tol(fm1), 1e-6)
  expect_equivalent(min_tol(fm2), 1e-2)
})
