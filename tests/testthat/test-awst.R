context("test-awst")

test_that("function gives consistent results", {
  old_retval <- matrix(data = c(-2, -2, 2.614586, -2, 2.541772,
                                2.395290, 2.253305, 2.391610, -2, -1.824604,
                                -2, -2, -1.999993, -1.867436, 2.596351,
                                -2, 2.662462, -1.999993, -2, -1.824604,
                                -2, -1.991209, -1.999993, -2, -1.824604),
                       nrow = 5, ncol = 5)

  set.seed(222)
  x <- matrix(rpois(25, lambda=5), ncol=5, nrow=5)
  retval <- awst(x)
  expect_true(max(abs(retval - old_retval)) < 1e-5)
})

