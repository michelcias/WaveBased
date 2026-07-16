# Manual benchmark: exact Daubechies-Lagarias evaluation vs table lookup with
# linear interpolation (wavelet.table argument).
#
# Run interactively; not executed by R CMD check.

library(WaveBased)

set.seed(123)

family <- "symmlets"
filter.size <- 20

cat("Building the interpolation table...\n")
print(system.time(tab <- wtable(family = family, filter.size = filter.size)))
print(tab)

results <- data.frame()

for(n in c(10^3, 10^4, 10^5)){
  x <- sort(runif(n))
  for(J in c(5, 7, 9)){
    t.exact <- system.time(
      PHI(x, J = J, family = family, filter.size = filter.size)
    )["elapsed"]
    t.interp <- system.time(
      PHI(x, J = J, family = family, filter.size = filter.size,
          wavelet.table = tab)
    )["elapsed"]
    results <- rbind(results,
                     data.frame(n = n, J = J, exact = t.exact,
                                interp = t.interp,
                                speedup = t.exact/max(t.interp, 1e-4)))
  }
}

rownames(results) <- NULL
cat("\nPHI: exact vs interpolated (seconds)\n")
print(results, digits = 3)

# Accuracy summary on the largest case
x <- sort(runif(10^5))
e <- max(abs(PHI(x, J = 9, family = family, filter.size = filter.size) -
             PHI(x, J = 9, family = family, filter.size = filter.size,
                 wavelet.table = tab)))
cat("\nMax abs difference at n = 1e5, J = 9:", format(e, digits = 3), "\n")
