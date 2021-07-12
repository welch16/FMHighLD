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

test_that("FMParam with 'all' strategy works", {

  fmparam <- FMParam(strategy = "all", prob = .9, k = 3)
  expect_equal(params_all(fmparam), list(prob = .9, k = 3))

})

test_that("FMParam with 'pick_m' strategy workds", {

  fmparam <- FMParam(strategy = "pick_m", prob = .9, m = 3, k = 3,
    m_sel = "random")
  expect_equal(params_pickm(fmparam),
    list(prob = 0.9, k = 3, m = 3, msel = "random"))

})