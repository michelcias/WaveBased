# Unit tests for the highest posterior density intervals.

test_that("hpdi returns the shortest window of the sorted sample", {

  set.seed(3)
  x <- rgamma(1000, shape = 2)

  got <- hpdi(x, prob = 0.9)
  s <- sort(x)
  m <- floor(0.9*1000)
  width <- s[(m + 1):1000] - s[1:(1000 - m)]
  best <- which.min(width)

  expect_equal(unname(got), c(s[best], s[best + m]))
  expect_named(got, c("lower", "upper"))

  # On a symmetric sample the interval is nearly centred, and its width is the
  # smallest one among the windows of the same span.
  z <- qnorm(seq(0.0005, 0.9995, length.out = 1000))
  eq <- hpdi(z, prob = 0.95)
  m <- floor(0.95*1000)
  expect_lt(abs(mean(eq)), 0.02)
  expect_equal(unname(diff(eq)), min(z[(m + 1):1000] - z[1:(1000 - m)]))
})


test_that("hpdi works column by column and validates its arguments", {

  set.seed(5)
  m <- matrix(rnorm(500*4), 500, 4)
  got <- hpdi(m, prob = 0.8)

  expect_equal(dim(got), c(4L, 2L))
  expect_equal(colnames(got), c("lower", "upper"))
  for(j in 1:4)
    expect_equal(got[j, ], hpdi(m[, j], prob = 0.8))

  expect_error(hpdi("a"), "numeric")
  expect_error(hpdi(rnorm(10), prob = 1), "between 0 and 1")
  expect_error(hpdi(1, prob = 0.9), "at least 2 draws")
})
