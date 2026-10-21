#==============================================================================
# Checks of a user-supplied filter (family = "Own")
#==============================================================================

# Checks that a user-supplied filter is the low-pass filter of an orthonormal
# wavelet basis of compact support, with the conditions that hold for any
# such filter: (i) it sums to sqrt(2); (ii) it is orthogonal to its
# translations by even shifts; (iii) Lawton's condition. (i) and (ii)
# together give H(pi) = 0, the one vanishing moment that the theory
# requires; how many more the filter has is the choice of its designer, and
# is not checked. (i) and (ii) are not enough: c(1, 0, 0, 1)/sqrt(2)
# satisfies both, but its scaling function is 1[0,3)/3, whose integer
# translations are not orthonormal. For a filter that satisfies (i) and
# (ii), the translations are orthonormal if, and only if, the eigenvalue 1
# of the transition matrix of the autocorrelation of the filter is simple
# (Lawton, 1991, J. Math. Phys. 32, 57-61).
#
# The coefficients may come with a few digits only, so (i) and (ii) are
# required up to 1e-6, with a warning above 1e-10 (the Symmlets of the
# package reach 1.4e-11). The eigenvalue 1 moves by about as much as the
# defect, hence the tolerance of the multiplicity. A second eigenvalue near 1
# means a filter close to failing Lawton's condition: its scaling function
# is very irregular and the algorithms that evaluate it converge slowly.
#
# Stops with the call of the exported function that checks the filter, and
# returns the filter as a double vector, invisibly.
.wb_check_filter <- function(wavelet.filter, call = sys.call(-1L)){
  h <- wavelet.filter
  if(!is.numeric(h) || length(h) < 2L || length(h) %% 2L != 0L ||
     any(!is.finite(h)))
    stop(simpleError("'wavelet.filter' must be a numeric vector of even length, with finite elements.",
                     call))
  h <- as.double(h)
  L <- length(h)

  d <- .wb_qmf_defect(h)
  if(d > 1e-6){
    flip <- (-1)^(seq_len(L) - 1L)*rev(h)
    hint <- if(.wb_qmf_defect(-h) <= 1e-6)
              " It sums to -sqrt(2): multiply it by -1."
            else if(.wb_qmf_defect(sqrt(2)*h) <= 1e-6)
              " It sums to 1: multiply it by sqrt(2)."
            else if(min(.wb_qmf_defect(flip), .wb_qmf_defect(-flip)) <= 1e-6)
              " It looks like the high-pass (wavelet) filter: provide the low-pass (scaling) one."
            else
              ""
    stop(simpleError(paste0("'wavelet.filter' is not the filter of an orthonormal wavelet basis: it must sum to sqrt(2) and be orthogonal to its translations by even shifts (largest defect ",
                            signif(d, 2), ").", hint), call))
  }
  if(d > 1e-10)
    warning(.wb_filter_warning(paste0("'wavelet.filter' satisfies the conditions of an orthonormal wavelet filter up to ",
                                      signif(d, 2), " only. The results carry an error of this order."),
                               call))

  ev <- eigen(.wb_lawton_matrix(h), only.values = TRUE)$values
  gap <- sort(abs(ev - 1))[2L]
  if(gap < max(1e-6, 100*d))
    stop(simpleError("'wavelet.filter' fails Lawton's condition: the integer translations of its scaling function are not orthonormal, although the filter sums to sqrt(2) and is orthogonal to its translations by even shifts (as c(1, 0, 0, 1)/sqrt(2) does).",
                     call))
  if(gap < 1e-3)
    warning(.wb_filter_warning(paste0("'wavelet.filter' nearly fails Lawton's condition (an eigenvalue at ",
                                      signif(gap, 2), " from 1): its scaling function is very irregular, and its evaluation may be inaccurate."),
                               call))

  invisible(h)
}

# Largest defect of the conditions on the coefficients of an orthonormal
# wavelet filter: |sum(h) - sqrt(2)| and |sum(h[n] h[n + 2k]) - delta_k|,
# k = 0, ..., L/2 - 1.
.wb_qmf_defect <- function(h){
  L <- length(h)
  d <- abs(sum(h) - sqrt(2))
  for(k in 0:(L %/% 2L - 1L)){
    i <- seq_len(L - 2L*k)
    d <- max(d, abs(sum(h[i]*h[i + 2L*k]) - (k == 0L)))
  }
  d
}

# Transition matrix of Lawton's condition, T[i, j] = a[2i - j] for i, j in
# -(L - 1), ..., L - 1, where a is the autocorrelation of the filter.
.wb_lawton_matrix <- function(h){
  L <- length(h)
  idx <- -(L - 1L):(L - 1L)
  a <- vapply(idx, function(m){
    i <- max(1L, 1L - m):min(L, L - m)
    sum(h[i]*h[i + m])
  }, numeric(1L))
  m <- outer(idx, idx, function(i, j) 2L*i - j)
  M <- matrix(0, 2L*L - 1L, 2L*L - 1L)
  ok <- abs(m) <= L - 1L
  M[ok] <- a[m[ok] + L]
  M
}

# The warnings of .wb_check_filter have their own class, so that a function
# that checked the filter at its entry can muffle the repetitions of its
# internal calls (.wb_quiet_filter).
.wb_filter_warning <- function(message, call)
  structure(class = c("wb_filter_warning", "warning", "condition"),
            list(message = message, call = call))

.wb_quiet_filter <- function(expr)
  withCallingHandlers(expr, wb_filter_warning = function(w)
    invokeRestart("muffleWarning"))
