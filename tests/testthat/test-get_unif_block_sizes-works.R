test_that("get_unif_block_sizes works", {

  nsnps <- 50
  nlds <- 11

  out <- rep(floor(nsnps / nlds), nlds)
  out[nlds] <- out[nlds] + nsnps - sum(out)

  expect_equal(get_unif_block_sizes(nsnps, nlds), out)

})
