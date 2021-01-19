
test_that("single trait formula works", {

  expect_equal(
    build_formula("z", c("a", "b"), FALSE),
    "z ~ 0 + a + b")

  expect_equal(
    build_formula("z", c("a", "b"), TRUE),
    "z ~ 1 + a + b")


})

test_that("multi trait formula works", {

  expect_equal(
    build_formula("z", c("a", "b"), FALSE,
      "c", "a", FALSE),
    "z ~ 0 + a + b + (0 + a || c)")

  expect_equal(
    build_formula("z", c("a", "b"), FALSE,
      "c", "a", TRUE),
    "z ~ 0 + a + b + (1 + a || c)")


})
