# ------------------------------------------------------------------------
# Generator of the precomputed Cohen-Daubechies-Vial boundary blocks
# (src/wav_filters_cdv.c).
#
# The blocks shipped with the package are derived here in multiple-precision
# arithmetic (Rmpfr, 320 bits), which sidesteps the ill conditioning of the
# boundary Gram matrices that limits the double-precision construction in
# R/cdv.R to moderate filter sizes. The resulting coefficients, rounded to
# double, define orthonormal bases accurate to machine precision for every
# tabulated filter.
#
# The derivation mirrors R/cdv.R step by step (same conventions, same block
# layout). The only structural difference is how the linear algebra is
# carried out in mpfr arithmetic:
#   * the restricted Gram matrix solves the Stein equation X = A X A' + C
#     by a doubling iteration instead of a large linear system;
#   * the edge scaling functions use an explicit Cholesky orthonormalization;
#   * the edge wavelets reuse the (well-conditioned) double-precision null
#     space structure from R/cdv.R and are then refined in mpfr by
#     projecting each vector onto the exact null space of its constraints.
#
# Each filter is validated, after rounding to double, by the same invariant
# used in the unit tests: the assembled one-step transform of the interval
# filter bank must be orthogonal to ~1e-14.
#
# Usage (from the package root):
#   Rscript --vanilla tools/make_cdv_tables.R
#
# Runtime is a few minutes; the script overwrites src/wav_filters_cdv.c.
# ------------------------------------------------------------------------

suppressMessages(library(Rmpfr))

PREC <- 320L                      # working precision in bits (~96 digits)
SIZES <- seq(4L, 24L, by = 2L)    # tabulated filter sizes per family
FAMILIES <- c(daublets = 1L, symmlets = 2L)

source("R/cdv.R")                 # double-precision reference machinery

# --- small mpfr linear-algebra helpers ------------------------------------

mzero <- function(nr, nc) mpfr2array(mpfr(rep(0, nr * nc), PREC), dim = c(nr, nc))

mnum <- function(x) mpfr(x, PREC)

# Gaussian elimination with partial pivoting, X = A^{-1} B.
msolve <- function(A, B) {
  n <- nrow(A)
  if (is.null(dim(B))) B <- mpfr2array(B, dim = c(n, 1L))
  m <- ncol(B)
  M <- mzero(n, n + m)
  M[, seq_len(n)] <- A
  M[, n + seq_len(m)] <- B
  for (k in seq_len(n)) {
    piv <- k - 1L + which.max(asNumeric(abs(M[k:n, k])))
    if (piv != k) {
      tmp <- M[k, ]; M[k, ] <- M[piv, ]; M[piv, ] <- tmp
    }
    M[k, ] <- M[k, ] / M[k, k]
    for (i in seq_len(n)[-k])
      if (asNumeric(M[i, k]) != 0)
        M[i, ] <- M[i, ] - M[i, k] * M[k, ]
  }
  M[, n + seq_len(m), drop = FALSE]
}

# Lower Cholesky factor of a symmetric positive definite mpfr matrix.
mchol_lower <- function(G) {
  n <- nrow(G)
  L <- mzero(n, n)
  for (j in seq_len(n)) {
    s <- G[j, j]
    if (j > 1) s <- s - sum(L[j, seq_len(j - 1)]^2)
    L[j, j] <- sqrt(s)
    if (j < n)
      for (i in (j + 1):n) {
        s <- G[i, j]
        if (j > 1) s <- s - sum(L[i, seq_len(j - 1)] * L[j, seq_len(j - 1)])
        L[i, j] <- s / L[j, j]
      }
  }
  L
}

# X = L^{-1} B for lower-triangular L.
mforwardsolve <- function(L, B) {
  n <- nrow(L)
  X <- mzero(n, ncol(B))
  for (i in seq_len(n)) {
    acc <- B[i, ]
    if (i > 1)
      for (j in seq_len(i - 1))
        acc <- acc - L[i, j] * X[j, ]
    X[i, ] <- acc / L[i, i]
  }
  X
}

# --- mpfr port of the derivation in R/cdv.R -------------------------------

# Orthogonality defect of the filter itself: the shipped coefficients are
# not all extended-precision (the Symmlets carry ~17 significant digits),
# and no construction can be more orthonormal than its filter. All the
# validation gates below adapt to this measured defect.
mqmf_defect <- function(h) {
  L <- length(h)
  d <- abs(sum(h) - sqrt(mnum(2)))
  for (k in 0:(L %/% 2 - 1)) {
    s <- sum(h[seq_len(L - 2 * k)] * h[seq_len(L - 2 * k) + 2 * k])
    d <- max(d, abs(s - (k == 0)))
  }
  asNumeric(d)
}

# Filter value helper on an mpfr filter vector (0-based index i).
mhval <- function(h, i) {
  out <- mpfr(rep(0, length(i)), PREC)
  ok <- i >= 0 & i <= length(h) - 1
  if (any(ok)) out[ok] <- h[i[ok] + 1]
  out
}

# Values of phi at the integers 1, ..., L-2: the eigenvector equation
# (sqrt(2) A - I) v = 0 is solved with the partition-of-unity normalization
# sum(v) = 1 replacing the last row.
mphi_integers <- function(h) {
  L <- length(h)
  ks <- seq_len(L - 2)
  A <- mzero(L - 2, L - 2)
  for (i in ks)
    for (j in ks)
      A[i, j] <- sqrt(mnum(2)) * mhval(h, 2L * i - j)
  M <- A
  for (i in ks) M[i, i] <- M[i, i] - 1
  M[L - 2, ] <- mpfr(rep(1, L - 2), PREC)
  rhs <- mzero(L - 2, 1L)
  rhs[L - 2, 1] <- mnum(1)
  as.vector(msolve(M, rhs))
}

# Moments of phi up to order nmom - 1 (M_0 = 1).
mmoments <- function(h, nmom) {
  L <- length(h)
  nn <- mpfr(0:(L - 1), PREC)
  M <- mpfr(rep(0, nmom), PREC)
  M[1] <- 1
  if (nmom > 1)
    for (beta in 1:(nmom - 1)) {
      acc <- mnum(0)
      for (gam in 0:(beta - 1))
        acc <- acc + choose(beta, gam) * sum(h * nn^(beta - gam)) * M[gam + 1]
      M[beta + 1] <- sqrt(mnum(2)) * 2^mnum(-beta - 1) * acc / (1 - 2^mnum(-beta))
    }
  M
}

# Restricted Gram matrix on the grid m = -(L-2), ..., 0. The straddling
# block solves the Stein equation X = Hs X Hs' + C by doubling:
# X = sum_k Hs^k C Hs'^k, with rho(Hs) = 2^(-1/2).
mrestricted_gram <- function(h) {
  L <- length(h)
  S <- -(L - 2):-1
  ns <- length(S)
  Hs <- mzero(ns, ns)
  Cc <- mzero(ns, ns)
  for (ip in seq_len(ns))
    for (ia in seq_len(ns))
      Hs[ip, ia] <- mhval(h, S[ia] - 2L * S[ip])
  for (ip in seq_len(ns))
    for (iq in seq_len(ns)) {
      arange <- (2L * S[ip]):(2L * S[ip] + L - 1L)
      brange <- (2L * S[iq]):(2L * S[iq] + L - 1L)
      common <- intersect(arange[arange >= 0], brange[brange >= 0])
      if (length(common))
        Cc[ip, iq] <- sum(mhval(h, common - 2L * S[ip]) *
                            mhval(h, common - 2L * S[iq]))
    }
  X <- Cc
  P <- Hs
  for (dbl in 1:10) {               # 2^10 terms: error ~ 2^(-1024)
    X <- X + P %*% X %*% t(P)
    P <- P %*% P
  }
  I0 <- mzero(ns + 1, ns + 1)
  I0[seq_len(ns), seq_len(ns)] <- X
  I0[ns + 1, ns + 1] <- mnum(1)
  I0
}

# Monomial coefficients of the shifted Legendre polynomials on [0, xmax].
mlegendre_shifted <- function(deg, xmax) {
  A <- mzero(deg + 1, deg + 1)
  A[1, 1] <- mnum(1)
  if (deg >= 1) {
    A[2, 2] <- mnum(1)
    if (deg >= 2)
      for (n in 1:(deg - 1))
        A[n + 2, ] <- ((2 * n + 1) * c(mnum(0), A[n + 1, seq_len(deg)]) -
                         n * A[n, ]) / (n + 1)
  }
  S <- mzero(deg + 1, deg + 1)
  for (r in 0:deg)
    for (s in 0:r)
      S[r + 1, s + 1] <- choose(r, s) * (2 / mnum(xmax))^s * (-1)^(r - s)
  A %*% S
}

# Edge scaling functions: Cholesky orthonormalization of the polynomial
# tails (in the shifted-Legendre basis) w.r.t. the restricted Gram matrix.
medge_scaling <- function(h, Igram) {
  L <- length(h)
  Nv <- L %/% 2
  mgrid <- mpfr(-(L - 2):0, PREC)
  Xmom <- mzero(Nv, L - 1)
  M <- mmoments(h, Nv)
  for (beta in 0:(Nv - 1)) {
    cf <- mpfr(rep(0, L - 1), PREC)
    for (gam in 0:beta)
      cf <- cf + choose(beta, gam) * mgrid^(beta - gam) * M[gam + 1]
    Xmom[beta + 1, ] <- cf
  }
  Cmat <- mlegendre_shifted(Nv - 1, L - 1) %*% Xmom
  G <- Cmat %*% Igram %*% t(Cmat)
  D <- mforwardsolve(mchol_lower(G), Cmat)
  err <- max(abs(asNumeric(D %*% Igram %*% t(D)) - diag(Nv)))
  stopifnot(err < 1e-50)
  D
}

# Two-scale blocks of the edge scaling functions (as in .cdv_twoscale).
mtwoscale <- function(h, D, Igram, tol) {
  L <- length(h)
  Nv <- L %/% 2
  mgrid <- -(L - 2):0
  ng <- length(mgrid)
  K <- mzero(ng, ng)
  for (ip in seq_len(ng))
    for (ia in seq_len(ng)) {
      ha <- mhval(h, mgrid[ia] - 2L * mgrid[ip])
      if (asNumeric(ha) != 0)
        K[ip, ] <- K[ip, ] + ha * Igram[ia, ]
    }
  B <- D %*% K %*% t(D)
  b <- mzero(Nv, L - 1)
  for (m in 1:(L - 1))
    b[, m] <- D %*% mhval(h, m - 2L * mgrid)
  err <- max(abs(asNumeric(B %*% t(B) + b %*% t(b)) - diag(Nv)))
  stopifnot(err < tol)
  list(B = B, b = b)
}

# Edge wavelets: the (well-conditioned) support structure and starting
# vectors come from the double-precision construction in R/cdv.R, applied
# to the rounded two-scale blocks; each vector is then refined in mpfr by
# projection onto the exact null space of its constraint matrix, keeping
# previously refined wavelets among the constraints so that the refined
# family is orthonormal to working precision.
medge_wavelets <- function(h, B, b, Bd, bd, tol) {
  L <- length(h)
  Nv <- L %/% 2
  hd <- asNumeric(h)
  dbl <- .cdv_edge_wavelets(hd, Bd, bd)
  uw <- max(ncol(dbl$u), L - 1)
  if (ncol(dbl$u) < uw)
    dbl$u <- cbind(dbl$u, matrix(0, Nv, uw - ncol(dbl$u)))
  g <- (-1)^(0:(L - 1)) * rev(h)
  P <- Nv + uw

  Cons0 <- mzero(Nv, P)             # coarse edge scaling functions
  Cons0[, seq_len(Nv)] <- B
  Cons0[, Nv + seq_len(L - 1)] <- b
  mp <- seq_len(uw %/% 2)
  Ch <- mzero(length(mp), P)        # overlapping interior scalings
  Cg <- mzero(length(mp), P)        # overlapping interior wavelets
  for (m1 in mp) {
    Ch[m1, Nv + seq_len(uw)] <- mhval(h, seq_len(uw) - 2L * m1)
    Cg[m1, Nv + seq_len(uw)] <- mhval(g, seq_len(uw) - 2L * m1)
  }

  U <- mzero(Nv, Nv)
  u <- mzero(Nv, uw)
  eps <- 2^mnum(-2 * PREC %/% 3)
  for (k in seq_len(Nv)) {
    rows <- rbind(Cons0, Ch, Cg)
    if (k > 1) {
      prev <- mzero(k - 1, P)
      prev[, seq_len(Nv)] <- U[seq_len(k - 1), , drop = FALSE]
      prev[, Nv + seq_len(uw)] <- u[seq_len(k - 1), , drop = FALSE]
      rows <- rbind(rows, prev)
    }
    K <- t(rows) %*% rows
    for (i in seq_len(P)) K[i, i] <- K[i, i] + eps
    w0 <- mpfr(c(dbl$U[k, ], dbl$u[k, ]), PREC)
    w <- as.vector(msolve(K, mpfr2array(w0, dim = c(P, 1L))))
    w <- w / sqrt(sum(w^2))
    if (asNumeric(sum(w * w0)) < 0) w <- -w
    U[k, ] <- w[seq_len(Nv)]
    u[k, ] <- w[Nv + seq_len(uw)]
  }
  err <- max(abs(asNumeric(U %*% t(U) + u %*% t(u)) - diag(Nv)))
  stopifnot(err < tol)
  list(U = U, u = u)
}

# Full one-side derivation in mpfr, returned already rounded to double.
mone_side <- function(h) {
  L <- length(h)
  tol <- max(1e-45, 1e3 * mqmf_defect(h))
  Igram <- mrestricted_gram(h)
  D <- medge_scaling(h, Igram)
  ts <- mtwoscale(h, D, Igram, tol)
  ew <- medge_wavelets(h, ts$B, ts$b, asNumeric(ts$B), asNumeric(ts$b), tol)
  phiint <- mphi_integers(h)
  pvec <- mpfr2array(c(rev(phiint), mnum(0)), dim = c(L - 1, 1L))
  phi0 <- asNumeric(D %*% pvec)
  list(B = asNumeric(ts$B), b = asNumeric(ts$b),
       U = asNumeric(ew$U), u = asNumeric(ew$u), phi0 = c(phi0))
}

# --- filter parsing (full precision, straight from the C sources) ---------

parse_filters_mpfr <- function(path) {
  txt <- readLines(path)
  filters <- list()
  cur <- NULL
  for (ln in txt) {
    m <- regmatches(ln, regexec("rfs == ([0-9]+)", ln))[[1]]
    if (length(m) == 2) cur <- m[2]
    m2 <- regmatches(ln, regexec("rwfilter\\[([0-9]+)\\][ ]*=[ ]*([-0-9.e+]+);", ln))[[1]]
    if (length(m2) == 3 && !is.null(cur)) {
      if (is.null(filters[[cur]])) filters[[cur]] <- character(0)
      filters[[cur]][as.integer(m2[2]) + 1] <- m2[3]
    }
  }
  lapply(filters, function(s) mpfr(s, PREC))
}

# --- generation and validation --------------------------------------------

derive_filter <- function(h) {
  left <- mone_side(h)
  right <- mone_side(rev(h))
  uwidth <- max(ncol(left$u), ncol(right$u))
  pad <- function(u) cbind(u, matrix(0, nrow(u), uwidth - ncol(u)))
  blocks <- list(BL = left$B, bL = left$b, UL = left$U, uL = pad(left$u),
                 phi0L = left$phi0,
                 BR = right$B, bR = right$b, UR = right$U, uR = pad(right$u),
                 phi0R = right$phi0, uwidth = as.integer(uwidth))

  # validate the rounded blocks with the unit-test invariant
  hd <- asNumeric(h)
  jmin <- .cdv_min_level(length(hd), uwidth, "decompose")
  O <- .cdv_step_matrix(hd, blocks, 2^(jmin + 1))
  err <- max(abs(O %*% t(O) - diag(2^(jmin + 1))))
  stopifnot(err < max(5e-14, 1e2 * mqmf_defect(h)))
  attr(blocks, "step.error") <- err
  blocks
}

fmt <- function(x) {
  out <- formatC(x, format = "e", digits = 17, width = 25)
  out[x == 0] <- formatC("0.0", width = 25)
  out
}

emit_side <- function(con, blocks, side) {
  vals <- c(blocks[[paste0("B", side)]], blocks[[paste0("b", side)]],
            blocks[[paste0("U", side)]], blocks[[paste0("u", side)]],
            blocks[[paste0("phi0", side)]])
  txt <- fmt(vals)
  for (i in seq(1, length(txt), by = 3))
    writeLines(paste0("  ", paste(txt[i:min(i + 2, length(txt))],
                                  collapse = ", "),
                      if (i + 2 < length(txt)) "," else ""), con)
}

main <- function() {
  filters <- list(
    daublets = parse_filters_mpfr("src/wav_filters_daublets.c"),
    symmlets = parse_filters_mpfr("src/wav_filters_symmlets.c"))

  con <- file("src/wav_filters_cdv.c", "w")
  on.exit(close(con))
  writeLines(c(
    "/**",
    " * @file wav_filters_cdv.c",
    " * @brief Precomputed Cohen-Daubechies-Vial boundary block coefficients.",
    " * @author Michel H. Montoril",
    " * @date 2026",
    " *",
    " * @details Generated by tools/make_cdv_tables.R -- do not edit by hand.",
    " *          The blocks are derived in 320-bit arithmetic (Rmpfr) and",
    " *          rounded to double, so every tabulated filter yields an",
    " *          orthonormal interval basis accurate to machine precision.",
    " *          Each entry stores, for one filter, the left and right edge",
    " *          blocks concatenated as B (Nv x Nv), b (Nv x (L-1)),",
    " *          U (Nv x Nv), u (Nv x uwidth) and phi0 (Nv), all column-major,",
    " *          with Nv = L/2. See R/cdv.R for the construction itself.",
    " */",
    "",
    "#include \"wav_filters_cdv.h\"",
    ""), con)

  registry <- character(0)
  for (famname in names(FAMILIES)) {
    famcode <- FAMILIES[[famname]]
    for (L in SIZES) {
      h <- filters[[famname]][[as.character(L)]]
      stopifnot(!is.null(h), length(h) == L)
      t0 <- proc.time()[3]
      blocks <- derive_filter(h)
      tag <- paste0(substr(famname, 1, 1), L)
      cat(sprintf("%-4s: uwidth=%2d  step.err=%8.1e  (%.1fs)\n", tag,
                  blocks$uwidth, attr(blocks, "step.error"),
                  proc.time()[3] - t0))
      writeLines(sprintf("static const double cdv_%s[] = {", tag), con)
      emit_side(con, blocks, "L")
      writeLines("  ,", con)
      emit_side(con, blocks, "R")
      writeLines("};", con)
      writeLines("", con)
      registry <- c(registry,
                    sprintf("  {%d, %d, %d, cdv_%s},", famcode, L,
                            blocks$uwidth, tag))
    }
  }

  writeLines("const cdv_table_entry cdv_table[] = {", con)
  writeLines(registry, con)
  writeLines("};", con)
  writeLines("", con)
  writeLines(sprintf("const int cdv_table_size = %d;", length(registry)), con)
  invisible(NULL)
}

main()
