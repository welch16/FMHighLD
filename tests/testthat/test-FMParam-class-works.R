test_that("FMParam default works", {
  expect_true(is(FMParam(), "FMParam"))
})

test_that("FMParam works when strategy == 'all'", {

  fmp <- FMParam(strategy = "all")
  expect_equal(strategy(fmp), "all")

})

test_that("FMParam works when strategy == 'pick_m'", {

  fmp <- FMParam(strategy = "pick_m")
  expect_equal(strategy(fmp), "pick_m")

})