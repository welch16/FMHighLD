test_that("FMParam default works", {

  fmparam <- FMParam()
  expect_equal(error_bound(fmparam), 1e-4)
  expect_equal(max_iter(fmparam), 100)
  expect_equal(min_tol(fmparam), 1e-6)
  expect_equal(annot_tol(fmparam), 1e-6)
  expect_equal(strategy(fmparam), "none")
})

test_that("FMParam with 'none' strategy works", {

  fmparam <- FMParam(strategy = "none")
  expect_equal(strategy(fmparam), "none")

})

