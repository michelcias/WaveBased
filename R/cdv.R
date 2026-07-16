#' @title Numerical Construction of the Cohen-Daubechies-Vial Boundary Blocks
#'
#' @description Internal machinery that derives, from an orthonormal wavelet
#' filter, all the coefficient blocks required by the boundary-corrected
#' (\emph{wavelets on the interval}) basis of Cohen, Daubechies and Vial (1993).
#' The construction is carried out numerically from first principles, so no
#' hard-coded coefficient tables are needed: everything follows from the
#' filter itself and is validated by internal orthonormality checks.
#'
#' The convention adopted here matches the rest of the package: the scaling
#' function \eqn{\phi} associated to a filter of length \eqn{L} has support
#' \eqn{[0, L-1]} and the wavelet filter has \eqn{N = L/2} vanishing moments.
#' On \eqn{[0, \infty)} the boundary-corrected multiresolution space keeps the
#' interior translates \eqn{\phi(\cdot - m)}, \eqn{m \ge 1}, and replaces the
#' translates with \eqn{m \le 0} by the \eqn{N} orthonormalized "polynomial
#' tails"
#' \deqn{T_\alpha = \sum_{m \le 0} c_m^\alpha \phi(\cdot - m)|_{[0,\infty)},
#' \quad c_m^\alpha = \langle x^\alpha, \phi(\cdot - m)\rangle,}
#' for \eqn{\alpha = 0, \ldots, N-1}, which is exactly the edge space of the
#' CDV construction (the resulting edge functions agree with the ones
#' published by CDV up to an orthogonal rotation within the edge block, which
#' leaves the multiresolution spaces, and hence all projections, unchanged).
#'
#' @references
#' Cohen, A., Daubechies, I. and Vial, P. (1993). Wavelets on the Interval and
#' Fast Wavelet Transforms. \emph{Applied and Computational Harmonic
#' Analysis}, 1(1), 54--81.
#'
#' @name cdv-internal
#' @keywords internal
NULL

# Environment used to cache the boundary blocks per filter, so the (cheap,
# but not free) linear algebra below runs only once per filter and session.
.cdv_cache <- new.env(parent = emptyenv())

# Values of the scaling function at the integers 1, ..., L-2, obtained as the
# eigenvector of the refinement matrix for eigenvalue 1, normalized so that
# sum(phi(k)) = 1 (partition of unity).
.cdv_phi_integers <- function(h) {
  L <- length(h)
  ks <- seq_len(L - 2)
  A <- outer(ks, ks, function(i, j) {
    idx <- 2 * i - j + 1
    val <- numeric(length(idx))
    ok <- idx >= 1 & idx <= L
    val[ok] <- h[idx[ok]]
    val
  })
  A <- sqrt(2) * A
  ei <- eigen(A)
  pos <- which.min(abs(ei$values - 1))
  if (abs(ei$values[pos] - 1) > 1e-8)
    stop("CDV construction: refinement matrix has no eigenvalue 1; ",
         "the provided filter does not look like an orthonormal wavelet filter.")
  v <- Re(ei$vectors[, pos])
  v / sum(v)
}

# Moments M[beta + 1] = integral of x^beta * phi(x) dx, beta = 0, ..., nmom-1,
# from the standard refinement recursion (with the normalization M_0 = 1).
.cdv_moments <- function(h, nmom) {
  L <- length(h)
  n <- 0:(L - 1)
  M <- numeric(nmom)
  M[1] <- 1
  if (nmom > 1) {
    for (beta in 1:(nmom - 1)) {
      acc <- 0
      for (gam in 0:(beta - 1))
        acc <- acc + choose(beta, gam) * sum(h * n^(beta - gam)) * M[gam + 1]
      M[beta + 1] <- sqrt(2) * 2^(-beta - 1) * acc / (1 - 2^(-beta))
    }
  }
  M
}

# Restricted Gram matrix I[p, q] = <phi(. - p), phi(. - q)> restricted to
# [0, Inf), for p, q = -(L-2), ..., 0. The unknown block (p, q <= -1) solves
# the linear system implied by the refinement equation; the row/column p = 0
# is exact (phi(. - 0) is entirely supported on [0, Inf)).
.cdv_restricted_gram <- function(h) {
  L <- length(h)
  S <- -(L - 2):-1          # straddling translates
  ns <- length(S)           # = L - 2
  hval <- function(i) ifelse(i >= 0 & i <= L - 1, h[pmin(pmax(i, 0), L - 1) + 1], 0)

  # System: I[p,q] = sum_{a,b} h[a-2p] h[b-2q] J[a,b], with
  # J[a,b] = I[a,b] if a,b in S; delta[a,b] if a,b >= 0; 0 otherwise.
  Tmat <- matrix(0, ns * ns, ns * ns)
  cvec <- numeric(ns * ns)
  amax <- L - 1              # h index range 0..L-1 => a in [2p, 2p+L-1]
  for (iq in seq_len(ns)) {
    q <- S[iq]
    for (ip in seq_len(ns)) {
      p <- S[ip]
      row <- ip + (iq - 1) * ns
      # constant part: sum over a = b >= 0
      arange <- seq.int(2 * p, 2 * p + amax)
      brange <- seq.int(2 * q, 2 * q + amax)
      common <- intersect(arange[arange >= 0], brange[brange >= 0])
      if (length(common))
        cvec[row] <- sum(hval(common - 2 * p) * hval(common - 2 * q))
      # unknown part: a, b in S
      for (ib in seq_len(ns)) {
        bb <- S[ib]
        hb <- hval(bb - 2 * q)
        if (hb == 0) next
        for (ia in seq_len(ns)) {
          aa <- S[ia]
          ha <- hval(aa - 2 * p)
          if (ha == 0) next
          Tmat[row, ia + (ib - 1) * ns] <- Tmat[row, ia + (ib - 1) * ns] + ha * hb
        }
      }
    }
  }
  x <- solve(diag(ns * ns) - Tmat, cvec)
  Istr <- matrix(x, ns, ns)
  Istr <- (Istr + t(Istr)) / 2   # enforce symmetry

  # Extend with the exact row/column for p = 0.
  I0 <- matrix(0, ns + 1, ns + 1)
  I0[seq_len(ns), seq_len(ns)] <- Istr
  I0[ns + 1, ns + 1] <- 1
  rownames(I0) <- colnames(I0) <- as.character(c(S, 0))
  I0
}

# Monomial coefficients of the Legendre polynomials of degree 0, ..., deg,
# shifted to the interval [0, xmax]. Row alpha + 1 holds the coefficients of
# the polynomial in increasing powers of x. Using these (instead of raw
# monomials) as the polynomial basis of the edge space keeps the edge Gram
# matrix well conditioned for longer filters.
.cdv_legendre_shifted <- function(deg, xmax) {
  A <- matrix(0, deg + 1, deg + 1)   # coefficients in t on [-1, 1]
  A[1, 1] <- 1
  if (deg >= 1) {
    A[2, 2] <- 1
    if (deg >= 2) {
      for (n in 1:(deg - 1)) {
        A[n + 2, ] <- ((2 * n + 1) * c(0, A[n + 1, seq_len(deg)]) -
                         n * A[n, ]) / (n + 1)
      }
    }
  }
  # substitute t = (2/xmax) x - 1
  S <- matrix(0, deg + 1, deg + 1)   # S[r+1, s+1] = coef of x^s in t^r
  for (r in 0:deg)
    for (s in 0:r)
      S[r + 1, s + 1] <- choose(r, s) * (2 / xmax)^s * (-1)^(r - s)
  A %*% S
}

# Orthonormalized edge scaling functions (left edge). Returns the Nv x (L-1)
# coefficient matrix D over the translate grid m = -(L-2), ..., 0, i.e.,
# phiL_k = sum_m D[k, m] phi(. - m) restricted to [0, Inf), obtained by a
# Gram-Schmidt (Cholesky) orthonormalization of the polynomial tails, in the
# order of increasing polynomial degree.
.cdv_edge_scaling <- function(h, Igram) {
  L <- length(h)
  Nv <- L %/% 2
  mgrid <- -(L - 2):0
  M <- .cdv_moments(h, Nv)
  # <x^beta, phi(. - m)> for beta = 0, ..., Nv - 1 on the translate grid
  Xmom <- matrix(0, Nv, length(mgrid))
  for (beta in 0:(Nv - 1)) {
    cf <- numeric(length(mgrid))
    for (gam in 0:beta)
      cf <- cf + choose(beta, gam) * mgrid^(beta - gam) * M[gam + 1]
    Xmom[beta + 1, ] <- cf
  }
  # tails of a well-conditioned polynomial basis (shifted Legendre on the
  # support [0, L-1] of the edge functions)
  Leg <- .cdv_legendre_shifted(Nv - 1, L - 1)
  Cmat <- Leg %*% Xmom
  # Iterated Cholesky orthonormalization w.r.t. the restricted Gram. The rows
  # of Cmat are exact, so every pass stays exactly inside the edge space (the
  # span is preserved); the iteration only improves orthonormality, which is
  # limited by the conditioning of the Gram matrix in double precision. The
  # final residual is the real quality gate.
  D <- Cmat
  err <- Inf
  for (pass in 1:4) {
    G <- D %*% Igram %*% t(D)
    G <- (G + t(G)) / 2
    R <- tryCatch(chol(G), error = function(e) NULL)
    if (is.null(R)) break
    D <- forwardsolve(t(R), D)
    err <- max(abs(D %*% Igram %*% t(D) - diag(Nv)))
    if (err < 1e-12) break
  }
  if (!is.finite(err) || err > 1e-9)
    stop("CDV construction: orthonormalization of the edge scaling functions ",
         "did not converge in double precision for this filter (residual ",
         format(err, digits = 3), "). Use a smaller filter.size with ",
         "boundary = \"interval\".")
  D
}

# Two-scale blocks of the left edge scaling functions:
#   phiL_k = sqrt(2) [ sum_l B[k,l] phiL_l(2.) + sum_{m=1}^{L-1} b[k,m] phi(2. - m) ].
# B is Nv x Nv and b is Nv x (L-1).
.cdv_twoscale <- function(h, D, Igram) {
  L <- length(h)
  Nv <- L %/% 2
  mgrid <- -(L - 2):0
  ng <- length(mgrid)
  hval <- function(i) ifelse(i >= 0 & i <= L - 1, h[pmin(pmax(i, 0), L - 1) + 1], 0)

  # K[p, q] = <phi(. - p)|[0,Inf), sqrt(2) phi(2. - q)|[0,Inf)>
  #         = sum_a h[a - 2p] I[a, q],  a restricted to the grid of Igram.
  K <- matrix(0, ng, ng)
  for (ip in seq_len(ng)) {
    p <- mgrid[ip]
    for (ia in seq_len(ng)) {
      a <- mgrid[ia]
      ha <- hval(a - 2 * p)
      if (ha == 0) next
      K[ip, ] <- K[ip, ] + ha * Igram[ia, ]
    }
  }
  B <- D %*% K %*% t(D)

  b <- matrix(0, Nv, L - 1)
  for (m in 1:(L - 1)) {
    hm <- hval(m - 2 * mgrid)
    b[, m] <- D %*% hm
  }

  err <- max(abs(B %*% t(B) + b %*% t(b) - diag(Nv)))
  if (err > 1e-8)
    stop("CDV construction: two-scale relation of the edge scaling functions ",
         "is not orthonormal (residual ", format(err, digits = 3), "). ",
         "The filter may be too long for a stable interval construction.")
  list(B = B, b = b)
}

# Left edge wavelets with (adaptively) staggered supports. The k-th edge
# wavelet lives, at the finer scale, on a local zone made of the Nv edge
# scaling functions plus the first M interior translates, for the smallest M
# admitting a new unit vector orthogonal to all coarse scaling functions, all
# overlapping interior wavelets and the previously constructed edge wavelets.
# Zones are grown one translate at a time, so filters whose minimal edge
# wavelet supports differ (they do across families) are handled uniformly.
# Returns U (Nv x Nv, edge part) and u (Nv x uwidth, interior part,
# zero-padded), where uwidth <= 2L - 3 is the largest zone actually used.
.cdv_edge_wavelets <- function(h, B, b) {
  L <- length(h)
  Nv <- L %/% 2
  g <- (-1)^(0:(L - 1)) * rev(h)   # g[n+1] = (-1)^n h[L-1-n]
  hval <- function(i) ifelse(i >= 0 & i <= L - 1, h[pmin(pmax(i, 0), L - 1) + 1], 0)
  gval <- function(i) ifelse(i >= 0 & i <= L - 1, g[pmin(pmax(i, 0), L - 1) + 1], 0)

  Mcap <- 4 * L                    # generous safety cap on the zone size
  U <- matrix(0, Nv, Nv)
  u <- matrix(0, Nv, Mcap)
  found <- 0

  for (M in 0:Mcap) {
    P <- Nv + M
    rows <- list(cbind(B, b[, seq_len(min(M, L - 1)), drop = FALSE],
                       matrix(0, Nv, max(0, M - (L - 1)))))
    mp <- seq_len(M %/% 2)
    if (length(mp)) {
      rows[[2]] <- t(vapply(mp, function(m1)
        c(rep(0, Nv), hval(seq_len(M) - 2 * m1)), numeric(P)))
      rows[[3]] <- t(vapply(mp, function(m1)
        c(rep(0, Nv), gval(seq_len(M) - 2 * m1)), numeric(P)))
    }
    if (found > 0)
      rows[[length(rows) + 1]] <- cbind(U[seq_len(found), , drop = FALSE],
                                        u[seq_len(found), seq_len(max(M, 1)),
                                          drop = FALSE][, seq_len(M),
                                                        drop = FALSE])
    Cons <- do.call(rbind, rows)

    sv <- svd(Cons, nu = 0, nv = P)
    dtol <- max(dim(Cons), P) * max(sv$d) * .Machine$double.eps * 100
    nullity <- P - sum(sv$d > dtol)
    if (nullity > 0) {
      W <- sv$v[, P - seq_len(nullity) + 1, drop = FALSE]
      for (jw in seq_len(ncol(W))) {
        if (found == Nv) break
        w <- W[, jw]
        w <- w * sign(w[which.max(abs(w))])   # deterministic sign
        found <- found + 1
        U[found, ] <- w[seq_len(Nv)]
        if (M > 0)
          u[found, seq_len(M)] <- w[Nv + seq_len(M)]
      }
    }
    if (found == Nv) {
      uwidth <- max(1L, max(which(colSums(abs(u)) > 1e-14), 0L))
      u <- u[, seq_len(uwidth), drop = FALSE]
      break
    }
  }
  if (found < Nv)
    stop("CDV construction: could not determine the edge wavelets for this ",
         "filter (found ", found, " of ", Nv, ").")

  err <- max(abs(U %*% t(U) + u %*% t(u) - diag(Nv)))
  if (err > 1e-8)
    stop("CDV construction: edge wavelets are not orthonormal (residual ",
         format(err, digits = 3), ").")
  list(U = U, u = u)
}

# Derive all blocks for one side from a filter (the right edge uses the
# reversed filter and index reflection, handled at evaluation/transform time).
.cdv_one_side <- function(h) {
  Igram <- .cdv_restricted_gram(h)
  D <- .cdv_edge_scaling(h, Igram)
  ts <- .cdv_twoscale(h, D, Igram)
  ew <- .cdv_edge_wavelets(h, ts$B, ts$b)
  phiint <- .cdv_phi_integers(h)          # phi(1), ..., phi(L-2)
  # phiL_k(0) = sum_p D[k, p] phi(-p), with -p = 0, ..., L-2 and phi(0) = 0.
  phi0 <- as.numeric(D %*% c(rev(phiint), 0))
  list(B = ts$B, b = ts$b, U = ew$U, u = ew$u, phi0 = phi0)
}

# Public (internal) entry point: all CDV blocks for a given filter, cached.
# The returned list layout is fixed and mirrored by the C code:
#   BL, bL, UL, uL, phi0L  -- left edge blocks (filter h)
#   BR, bR, UR, uR, phi0R  -- right edge blocks (reversed filter)
.cdv_blocks <- function(h) {
  key <- paste(format(h, digits = 17), collapse = ",")
  hit <- .cdv_cache[[key]]
  if (!is.null(hit))
    return(hit)
  if (length(h) %% 2 != 0 || length(h) < 4)
    stop("boundary = \"interval\" requires a filter of even length >= 4.")
  left <- .cdv_one_side(h)
  right <- .cdv_one_side(rev(h))
  # pad the interior wavelet blocks of the two sides to a common width
  uwidth <- max(ncol(left$u), ncol(right$u))
  pad <- function(u) cbind(u, matrix(0, nrow(u), uwidth - ncol(u)))
  out <- list(BL = left$B, bL = left$b, UL = left$U, uL = pad(left$u),
              phi0L = left$phi0,
              BR = right$B, bR = right$b, UR = right$U, uR = pad(right$u),
              phi0R = right$phi0,
              uwidth = as.integer(uwidth))
  assign(key, out, envir = .cdv_cache)
  out
}

# Minimum level for the interval basis. The level-j scaling basis needs the
# two edge zones not to interact (2^j >= 2L - 2, "basis"); the level-j
# wavelets are built inside level j+1 and need 2^(j+1) >= L + 2*uwidth so
# that the edge wavelet zones of the two sides stay disjoint ("wavelet");
# a decomposed basis with coarsest level j needs both ("decompose"). Since
# uwidth <= 2L - 3, the wavelet bound is at worst 2^(j+1) >= 5L - 6.
.cdv_min_level <- function(L, uwidth,
                           what = c("decompose", "basis", "wavelet")) {
  what <- match.arg(what)
  jbasis <- ceiling(log2(2 * L - 2))
  jwave <- ceiling(log2(L + 2 * uwidth)) - 1
  as.integer(switch(what,
                    basis = jbasis,
                    wavelet = jwave,
                    decompose = max(jbasis, jwave)))
}

# Filter coefficients of a built-in family (1 = Daublets, 2 = Symmlets,
# 3 = Coiflets), fetched from the C tables.
.wb_filter <- function(fam, filter.size) {
  .Call("_WaveBased_C_GetFilter", as.integer(fam), as.integer(filter.size))
}

# Resolves the boundary treatment requested through the 'boundary'/'periodic'
# arguments of the exported functions into the integer code used by the C
# entry points: 0 = raw ("none"), 1 = "periodic", 2 = "interval" (CDV).
.wb_boundary_code <- function(boundary, periodic, periodic_given) {
  if (is.null(boundary)) {
    boundary <- if (periodic) "periodic" else "none"
  } else {
    if (!is.character(boundary) || length(boundary) != 1L)
      stop("boundary must be one of \"periodic\", \"none\" or \"interval\".")
    boundary <- match.arg(tolower(boundary),
                          c("periodic", "none", "interval"))
    if (periodic_given && boundary != (if (periodic) "periodic" else "none"))
      warning("Argument 'boundary' overrides 'periodic'.")
  }
  match(boundary, c("none", "periodic", "interval")) - 1L
}

# Shared validation + block preparation for boundary = "interval". 'level'
# is the level whose minimum requirement must be met ('what' as in
# .cdv_min_level, with the same meaning of the level: j0 for a decomposed
# basis, J for PHI/PSI). Returns the CDV block list to be passed to C.
.cdv_prepare <- function(x, fam, filter.size, wavelet.filter, level, what) {
  if (fam == 3)
    stop("boundary = \"interval\" is not available for Coiflets: the ",
         "Cohen-Daubechies-Vial construction requires minimal-length ",
         "filters. Use Daublets or Symmlets.")
  h <- if (fam == 4) as.double(wavelet.filter)
       else .wb_filter(fam, filter.size)
  xf <- x[is.finite(x)]
  if (length(xf) && (min(xf) < 0 || max(xf) > 1))
    stop("boundary = \"interval\" requires the data to lie in [0, 1]. ",
         "Rescale x before calling this function.")
  blocks <- .cdv_blocks(h)
  jmin <- .cdv_min_level(length(h), blocks$uwidth, what)
  if (level < jmin)
    stop("boundary = \"interval\" with a filter of size ", length(h),
         " requires ", if (what == "decompose") "j0" else "J", " >= ", jmin,
         " (got ", level, "). Increase the level or reduce filter.size.")
  blocks
}

# Reference (pure R) one-step analysis matrix, mapping the 2^(j+1) basis
# values/coefficients to [scaling level j | detail level j]. Used by the unit
# tests to validate the blocks and the index mapping against the C kernel.
.cdv_step_matrix <- function(h, blocks, n_in) {
  L <- length(h)
  Nv <- L %/% 2
  uw <- ncol(blocks$uL)
  g <- (-1)^(0:(L - 1)) * rev(h)
  nout <- n_in %/% 2
  nint_in <- n_in - L            # interior translates m = 1, ..., nint_in
  nint_out <- nout - L
  stopifnot(nint_out >= 0, nint_in >= 2 * uw)
  O <- matrix(0, n_in, n_in)
  # index helpers on the input layout [edgeL | interior m | edgeR]
  icol <- function(m) Nv + m                    # interior m >= 1
  ircol <- function(mt) Nv + (nint_in + 1 - mt) # mirrored interior index
  # --- scaling rows ---
  for (k in seq_len(Nv)) {
    O[k, seq_len(Nv)] <- blocks$BL[k, ]
    O[k, icol(seq_len(L - 1))] <- blocks$bL[k, ]
  }
  for (m1 in seq_len(nint_out)) {
    O[Nv + m1, icol(2 * m1 + 0:(L - 1))] <- h
  }
  for (k in seq_len(Nv)) {
    O[nout - Nv + k, Nv + nint_in + seq_len(Nv)] <- blocks$BR[k, ]
    O[nout - Nv + k, ircol(seq_len(L - 1))] <- blocks$bR[k, ]
  }
  # --- detail rows ---
  for (k in seq_len(Nv)) {
    O[nout + k, seq_len(Nv)] <- blocks$UL[k, ]
    O[nout + k, icol(seq_len(uw))] <- blocks$uL[k, ]
  }
  for (m1 in seq_len(nint_out)) {
    O[nout + Nv + m1, icol(2 * m1 + 0:(L - 1))] <- g
  }
  for (k in seq_len(Nv)) {
    O[nout + nout - Nv + k, Nv + nint_in + seq_len(Nv)] <- blocks$UR[k, ]
    O[nout + nout - Nv + k, ircol(seq_len(uw))] <- blocks$uR[k, ]
  }
  O
}
