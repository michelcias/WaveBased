# Unit tests for the Polya-Gamma sampler that the logit link of bwregime() is
# built on.

rpg1 <- function(z) .Call("_WaveBased_C_rpg1", as.double(z))

# The first two moments of PG(1, z), which are known in closed form.
pg_mean <- function(z)
  ifelse(abs(z) < 1e-8, 0.25, tanh(z/2)/(2*z))

pg_var <- function(z)
  ifelse(abs(z) < 1e-8, 1/24, (sinh(z) - z)/(2*z^3*(1 + cosh(z))))

# The infinite sum representation, truncated. It is a second implementation of
# the same distribution, and an independent one: the exact sampler is a
# rejection scheme and this is a sum of exponential variables.
rpg1_sum <- function(z, m, K = 500){
  g <- matrix(stats::rexp(K*m), K, m)
  colSums(g/((1:K - 0.5)^2 + z^2/(4*pi^2)))/(2*pi^2)
}


test_that("the draws have the moments of the Polya-Gamma distribution", {

  set.seed(20)

  for(z in c(0, 0.5, 1, 2, 5, 10, 30)){

    x <- rpg1(rep(z, 40000))

    expect_true(all(x > 0))

    # Three standard errors of the sample mean, which is the tolerance a test
    # of a Monte Carlo quantity can afford without becoming flaky.
    expect_lt(abs(mean(x) - pg_mean(z)), 3*sqrt(pg_var(z)/40000))
    expect_lt(abs(var(x)/pg_var(z) - 1), 0.05)
  }
})


test_that("the exact sampler agrees with the infinite sum representation", {

  set.seed(21)

  for(z in c(0, 1, 5, 20)){

    x <- rpg1(rep(z, 4000))
    y <- rpg1_sum(z, 4000)

    expect_gt(suppressWarnings(stats::ks.test(x, y)$p.value), 0.01)
  }
})


test_that("the sampler is reproducible and takes the sign of its argument", {

  set.seed(22)
  a <- rpg1(c(-3, 0, 3, 100))

  set.seed(22)
  b <- rpg1(c(3, 0, -3, -100))

  # The distribution depends on the tilting parameter through its absolute
  # value, so the same stream gives the same draws.
  expect_equal(a, b)
  expect_true(all(is.finite(a)))

  # A large tilting parameter concentrates the variable near zero, and the
  # sampler must not underflow to it.
  expect_true(all(rpg1(rep(500, 100)) > 0))
})
