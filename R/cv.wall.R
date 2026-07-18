#' @title Cross-validation for the wall classifier
#'
#' @description The \code{cv.wall} function performs k-fold cross-validation
#' to jointly select the resolution level \code{J} and the penalty parameter
#' \eqn{\lambda} of the \command{\link{wall}} classifier over a grid. For
#' each candidate value of \code{J}, the whole path of \eqn{\lambda} is
#' cross-validated by \command{\link[glmnet]{cv.glmnet}}, using the same
#' folds; the pair minimizing the cross-validated measure is then selected.
#'
#' @param x Matrix (or data frame) of covariates, with \eqn{n} rows
#'   (observations) and \eqn{d} columns (variables). A vector is interpreted
#'   as a single covariate. For the \code{print} method, \code{x} is an
#'   object of class \code{"cv.wall"}.
#' @param y Response with two classes, as in \command{\link{wall}}.
#' @param J Vector with the grid of candidate resolution levels (single
#'   integers, used for all the covariates). By default, the grid
#'   \code{(j0+1):max(j0+2, ceiling(log2(n)/2))} is used.
#' @param j0 The coarsest resolution level of the decomposed wavelet basis.
#'   Default is \code{j0 = 0}. See \command{\link{wall}}.
#' @param family,filter.size,prec.wavelet,wavelet.filter The wavelet basis
#'   specification, as in \command{\link{wall}}.
#' @param boundary The boundary treatment of the basis, as in
#'   \command{\link{wall}}.
#' @param rescale Logical. If \code{TRUE} (default), the covariates are
#'   rescaled to \eqn{[}\code{eps}\eqn{, 1-}\code{eps}\eqn{]} before the
#'   basis evaluation. See \command{\link{wall}}.
#' @param eps A value in \eqn{[0, 0.5)} used by the rescaling. By default,
#'   \code{eps} \eqn{= 1.9^{-J}} for each candidate \code{J} (periodized
#'   basis); with \code{share.design = TRUE}, a single value
#'   \eqn{1.9^{-\max(J)}} is used for the whole grid. See
#'   \command{\link{wall}}.
#' @param type.measure The loss used for cross-validation, among the measures
#'   available in \command{\link[glmnet]{cv.glmnet}} for logistic regression:
#'   \code{"class"} (default, the misclassification error, as adopted in the
#'   numerical studies of Montoril, 2026), \code{"deviance"} (the binomial
#'   deviance) or \code{"auc"} (the area under the ROC curve, which is
#'   maximized instead of minimized).
#' @param nfolds Number of folds. Default is \code{nfolds = 10}. Ignored if
#'   \code{foldid} is provided.
#' @param foldid Optional vector with the fold membership (values in
#'   \code{1, ..., nfolds}) of each observation. By default, the folds are
#'   sampled once and shared by all the candidate values of \code{J}, so that
#'   the resolution levels are compared under the same data partitions.
#' @param use.table,wavelet.table Fast evaluation of the basis functions by
#'   table lookup, as in \command{\link{wall}}. With \code{use.table =
#'   "auto"} (default), the table is built once and reused across the whole
#'   grid of \code{J}, where it typically pays off.
#' @param sparse Sparse storage of the design matrix, as in
#'   \command{\link{wall}}.
#' @param share.design Logical. If \code{TRUE}, the design matrix of basis
#'   functions is built only once, at the largest candidate resolution
#'   level, and every candidate \code{J} reuses the corresponding subset of
#'   its columns -- the multiresolution spaces are nested, so the design of
#'   each \code{J} is exactly the leading block of columns of the design of
#'   \code{max(J)}, per covariate. This divides the cost of the basis
#'   evaluation by roughly the size of the grid. It requires a single
#'   rescaling constant \code{eps} for the whole grid (otherwise the
#'   rescaled covariates would change with \code{J}): if \code{eps} is not
#'   provided, \eqn{1.9^{-\max(J)}} is used. The default is \code{FALSE},
#'   which reproduces the tuning protocol of Montoril (2026) exactly
#'   (\code{eps} \eqn{= 1.9^{-J}} recomputed for each candidate).
#' @param lambda Optional user-supplied sequence of values of \eqn{\lambda},
#'   shared by all the candidate values of \code{J}. By default, a sequence
#'   is computed for each \code{J}; see \command{\link[glmnet]{glmnet}}.
#' @param standardize Logical. Should the columns of basis functions be
#'   internally standardized by \pkg{glmnet}? Default is \code{FALSE}. See
#'   \command{\link{wall}}.
#' @param weights Optional vector of observation weights.
#' @param parallel Logical. If \code{TRUE}, the cross-validation of each
#'   candidate \code{J} is parallelized over the folds; it requires a
#'   parallel backend registered for \pkg{foreach}, as in
#'   \command{\link[glmnet]{cv.glmnet}}. Ignored (with a warning) when
#'   \code{parallel.J = TRUE}.
#' @param parallel.J Logical. If \code{TRUE}, the grid of candidate
#'   resolution levels is processed in parallel, with one task per value of
#'   \code{J} -- each task builds its design matrix and cross-validates its
#'   whole path of \eqn{\lambda}. It requires a parallel backend registered
#'   for \pkg{foreach} (e.g., \pkg{doParallel}). The candidates are
#'   independent given the folds, so the results are identical to the
#'   sequential ones. This usually amortizes the work better than
#'   \code{parallel}, which parallelizes only the folds within each
#'   candidate. Default is \code{FALSE}.
#' @param trace Logical. If \code{TRUE}, prints a progress line for each
#'   candidate value of \code{J}.
#' @param ... Further arguments passed to
#'   \command{\link[glmnet]{cv.glmnet}} and, from there, to
#'   \command{\link[glmnet]{glmnet}}. For the \code{print} method, further
#'   arguments are ignored.
#'
#' @details
#' The procedure mirrors the tuning strategy in the numerical studies of
#' Montoril (2026): the folds are drawn once, each candidate resolution
#' level \code{J} has its full LASSO path cross-validated under those folds,
#' and the selected classifier is the one whose pair
#' \eqn{(J, \lambda)} attains the best cross-validated measure. Within the
#' selected \code{J}, both \code{lambda.min} (the minimizer) and
#' \code{lambda.1se} (the largest \eqn{\lambda} whose measure is within one
#' standard error of the minimum) are reported, as in
#' \command{\link[glmnet]{cv.glmnet}}.
#'
#' For each candidate \code{J}, the design matrix of basis functions is
#' built only once on the full data set (with the rescaling constant
#' \code{eps} \eqn{= 1.9^{-J}} recomputed accordingly) and reused across the
#' folds. When the basis is evaluated by table lookup (see
#' \command{\link{wall}}), a single table serves the whole grid of \code{J},
#' which makes the cross-validation substantially faster.
#'
#' With \code{share.design = TRUE}, the design is built only once for the
#' whole grid, at \code{max(J)}, and each candidate reuses the appropriate
#' subset of its columns; since the rescaling must then be the same for all
#' the candidates, a single \code{eps} (by default \eqn{1.9^{-\max(J)}}) is
#' used instead of \eqn{1.9^{-J}}. For a fixed common \code{eps}, the
#' subsetting is exact: the results are identical to the ones obtained by
#' rebuilding the design at each \code{J}.
#'
#' The final fit, stored in the component \code{wall.fit}, reuses the model
#' already fitted on the full data set during the cross-validation of the
#' selected \code{J} -- no refit is needed.
#'
#' @return An object of class \code{"cv.wall"}, i.e., a list with components
#'   \item{J}{The grid of candidate resolution levels.}
#'   \item{cv}{A list with one element per candidate \code{J}, each
#'   containing the vectors \code{lambda}, \code{cvm}, \code{cvsd},
#'   \code{cvup}, \code{cvlo} and \code{nzero} of the cross-validation, and
#'   the values \code{lambda.min}, \code{lambda.1se}, \code{cvm.min},
#'   \code{cvsd.min} and \code{nzero.min} attained at \code{lambda.min}.}
#'   \item{cvtab}{A data frame summarizing the grid: for each \code{J}, the
#'   value of \code{lambda.min}, the cross-validated measure, its standard
#'   error and the number of non-zero coefficients.}
#'   \item{J.min}{The selected resolution level.}
#'   \item{lambda.min, lambda.1se}{The selected values of \eqn{\lambda}
#'   within \code{J.min}.}
#'   \item{cvm.min, cvsd.min, nzero.min}{The cross-validated measure, its
#'   standard error and the number of non-zero coefficients at the selected
#'   pair.}
#'   \item{type.measure, measure.name}{The loss used.}
#'   \item{nfolds, foldid}{The folds used.}
#'   \item{wall.fit}{The selected classifier, an object of class
#'   \code{"wall"} fitted on the full data set with \code{J = J.min}. See
#'   \command{\link{wall}}.}
#'   \item{call}{The matched call.}
#'
#' @references
#' Montoril, M. H. (2026). Wavelet-based classification: an approach via
#' additive logistic regression. \emph{Manuscript}.
#'
#' Friedman, J., Hastie, T. and Tibshirani, R. (2010). Regularization Paths
#' for Generalized Linear Models via Coordinate Descent. \emph{Journal of
#' Statistical Software}, 33(1), 1--22, \doi{10.18637/jss.v033.i01}.
#'
#' @seealso \command{\link{wall}}, \command{\link{predict.cv.wall}},
#'   \command{\link{plot.cv.wall}}, \command{\link[glmnet]{cv.glmnet}}
#'
#' @author Michel H. Montoril \email{michel@@ufscar.br}
#'
#' @examples
#' # Simulated additive logistic model with two covariates
#' set.seed(123)
#' n <- 400
#' x <- matrix(runif(2*n), n, 2)
#' h <- 3*sin(2*pi*x[, 1]) + 8*(x[, 2] - 0.5)   # log-odds
#' y <- rbinom(n, 1, 1/(1 + exp(-h)))
#'
#' # Cross-validated choice of J and lambda
#' set.seed(1)
#' cvfit <- cv.wall(x, y, J = 1:4, filter.size = 8, nfolds = 5)
#' cvfit
#' plot(cvfit)
#'
#' # Predictions with the selected pair (J, lambda)
#' yhat <- predict(cvfit, x, type = "class")
#' mean(yhat != y)  # training misclassification rate
#'
#' @keywords classif
#' @importFrom glmnet cv.glmnet
#' @export
cv.wall <- function(x, y, J = NULL, j0 = 0, family = "Daublets",
                    filter.size = 20, prec.wavelet = 30, wavelet.filter,
                    boundary = c("periodic", "interval"),
                    rescale = TRUE, eps,
                    type.measure = c("class", "deviance", "auc"),
                    nfolds = 10, foldid = NULL,
                    use.table = c("auto", "always", "never"), wavelet.table = NULL,
                    sparse = c("auto", "always", "never"), share.design = FALSE,
                    lambda = NULL, standardize = FALSE, weights = NULL,
                    parallel = FALSE, parallel.J = FALSE, trace = FALSE, ...){

  this.call <- match.call()
  boundary <- match.arg(tolower(boundary[1L]), c("periodic", "interval"))
  use.table <- match.arg(use.table)
  sparse <- match.arg(sparse)
  type.measure <- match.arg(type.measure)

  x <- .wall_x(x)
  n <- nrow(x)
  d <- ncol(x)

  resp <- .wall_response(y)
  if(length(resp$y) != n)
    stop("The number of observations in 'y' must match the number of rows of 'x'.")

  j0 <- .wall_j0(j0)
  if(is.null(J))
    J <- seq(j0 + 1L, max(j0 + 2L, ceiling(log2(n)/2)))
  if(!is.numeric(J) || any(!is.finite(J)) || any(J != round(J)))
    stop("'J' must contain only integer values.")
  if(any(J <= j0))
    stop("'J' must be larger than 'j0' (= ", j0, ").")
  J <- sort(unique(as.integer(J)))

  if(missing(eps))
    eps <- NULL
  else if(length(eps) != 1L)
    stop("In cv.wall, 'eps' must be a single value; it is recycled over the covariates.")

  if(!is.logical(share.design) || length(share.design) != 1L || is.na(share.design))
    stop("'share.design' must be TRUE or FALSE.")
  if(!is.logical(parallel.J) || length(parallel.J) != 1L || is.na(parallel.J))
    stop("'parallel.J' must be TRUE or FALSE.")
  if(parallel.J && parallel){
    warning("'parallel' was ignored: with parallel.J = TRUE, the folds of each candidate J run sequentially inside its own task.")
    parallel <- FALSE
  }
  # A shared design requires a rescaling that does not change with J: fix
  # eps at its value for the largest candidate (Montoril et al., 2019 rule).
  if(share.design && is.null(eps) && rescale && boundary == "periodic")
    eps <- 1.9^(-max(J))

  wf.own <- if(missing(wavelet.filter)) NULL else wavelet.filter

  if(is.null(foldid))
    foldid <- sample(rep_len(seq_len(nfolds), n))
  else{
    if(length(foldid) != n)
      stop("'foldid' must have one entry per observation.")
    nfolds <- max(foldid)
  }
  if(nfolds < 3L)
    stop("'nfolds' must be at least 3.")

  # One lookup table serves the whole grid of J (and all the covariates).
  wtab <- .wall_table(use.table, wavelet.table, boundary,
                      workload = n*d*length(J),
                      family = family, filter.size = filter.size,
                      prec.wavelet = prec.wavelet, wavelet.filter = wf.own)

  drop.phi <- boundary == "periodic" && j0 == 0L
  L <- if(is.null(wf.own)) filter.size else length(wf.own)

  # The column ranges do not depend on J (only eps does), so they are
  # computed once and shared by the whole grid.
  xranges <- if(rescale) apply(x, 2L, range) else NULL

  # With share.design = TRUE, the design is built once at max(J); since eps
  # is fixed for the whole grid, the design of each candidate J is exactly
  # the leading 2^J (- drop.phi) columns of each covariate block.
  design.max <- NULL
  if(share.design){
    Jm <- rep_len(max(J), d)
    em <- .wall_eps(eps, Jm, j0, boundary, rescale)
    rsm <- .wall_rescale_pars(x, em, rescale, boundary, ranges = xranges)
    spec.max <- list(J = Jm, j0 = j0, boundary = boundary, family = family,
                     filter.size = filter.size, prec.wavelet = prec.wavelet,
                     wavelet.filter = wf.own, wavelet.table = wtab,
                     rescale = rescale, location = rsm$location,
                     scale = rsm$scale, eps = em, drop.phi = drop.phi,
                     sparse = .wall_sparse(sparse, Jm, j0, L, drop.phi),
                     xnames = colnames(x), classnames = resp$classnames,
                     nobs = n)
    design.max <- .wall_design(x, spec.max, clip = FALSE)
    blk.w <- 2^max(J) - as.integer(drop.phi)   # columns per covariate block
  }

  # All the work of one candidate J: design (built or sliced), whole-path
  # cross-validation, and the summaries needed afterwards. The candidates
  # are independent given foldid, so this closure can run sequentially
  # (lapply) or in parallel over the grid (foreach, parallel.J = TRUE)
  # with identical results.
  run.one <- function(i){

    Ji <- rep_len(J[i], d)
    ei <- .wall_eps(eps, Ji, j0, boundary, rescale)
    rs <- .wall_rescale_pars(x, ei, rescale, boundary, ranges = xranges)

    spec <- list(J = Ji, j0 = j0, boundary = boundary, family = family,
                 filter.size = filter.size, prec.wavelet = prec.wavelet,
                 wavelet.filter = wf.own, wavelet.table = wtab,
                 rescale = rescale, location = rs$location, scale = rs$scale,
                 eps = ei, drop.phi = drop.phi,
                 sparse = .wall_sparse(sparse, Ji, j0, L, drop.phi),
                 xnames = colnames(x), classnames = resp$classnames, nobs = n)

    if(share.design){
      wi <- 2^J[i] - as.integer(drop.phi)
      cols <- as.vector(outer(seq_len(wi), (seq_len(d) - 1L)*blk.w, `+`))
      design <- design.max[, cols, drop = FALSE]
    }
    else
      design <- .wall_design(x, spec, clip = FALSE)
    pf <- .wall_penalty(Ji, j0, drop.phi)

    cvfit <- glmnet::cv.glmnet(x = design, y = resp$y, family = "binomial",
                               weights = weights, lambda = lambda,
                               type.measure = type.measure, foldid = foldid,
                               penalty.factor = pf, standardize = standardize,
                               parallel = parallel, ...)

    imin <- match(cvfit$lambda.min, cvfit$lambda)
    res <- list(cv = list(J = J[i], lambda = cvfit$lambda, cvm = cvfit$cvm,
                          cvsd = cvfit$cvsd, cvup = cvfit$cvup,
                          cvlo = cvfit$cvlo, nzero = cvfit$nzero,
                          lambda.min = cvfit$lambda.min,
                          lambda.1se = cvfit$lambda.1se,
                          cvm.min = cvfit$cvm[imin],
                          cvsd.min = cvfit$cvsd[imin],
                          nzero.min = cvfit$nzero[imin]),
                name = cvfit$name, glmnet.fit = cvfit$glmnet.fit,
                spec = spec, pf = pf, nvars = ncol(design))

    if(trace && !parallel.J)
      cat(sprintf("J = %d: %s = %.5f at lambda = %.5g (%d nonzero)\n",
                  J[i], names(cvfit$name), res$cv$cvm.min,
                  res$cv$lambda.min, res$cv$nzero.min))

    res
  }

  if(parallel.J){
    if(!requireNamespace("foreach", quietly = TRUE))
      stop("'parallel.J = TRUE' requires the 'foreach' package.")
    `%dopar%` <- foreach::`%dopar%`
    runs <- foreach::foreach(i = seq_along(J)) %dopar% run.one(i)
  }
  else
    runs <- lapply(seq_along(J), run.one)

  cvlist <- lapply(runs, `[[`, "cv")

  if(trace && parallel.J)
    for(i in seq_along(J))
      cat(sprintf("J = %d: %s = %.5f at lambda = %.5g (%d nonzero)\n",
                  J[i], names(runs[[i]]$name), cvlist[[i]]$cvm.min,
                  cvlist[[i]]$lambda.min, cvlist[[i]]$nzero.min))

  # Selection of the best candidate (first minimum/maximum on ties, as in
  # the sequential comparison used before).
  cvm.all <- vapply(cvlist, `[[`, 0, "cvm.min")
  best.i <- if(type.measure == "auc") which.max(cvm.all) else which.min(cvm.all)
  best <- cvlist[[best.i]]
  best.name <- runs[[best.i]]$name

  # The full-data fit of the selected J, reused as the final classifier.
  best.spec <- runs[[best.i]]$spec
  best.spec$glmnet.fit <- runs[[best.i]]$glmnet.fit
  best.spec$lambda <- best.spec$glmnet.fit$lambda
  best.spec$penalty.factor <- runs[[best.i]]$pf
  best.spec$nvars <- runs[[best.i]]$nvars

  best.spec$call <- this.call
  class(best.spec) <- "wall"

  cvtab <- data.frame(J = J,
                      lambda.min = vapply(cvlist, `[[`, 0, "lambda.min"),
                      measure = vapply(cvlist, `[[`, 0, "cvm.min"),
                      sd = vapply(cvlist, `[[`, 0, "cvsd.min"),
                      nzero = vapply(cvlist, function(z) as.integer(z$nzero.min),
                                     0L))
  names(cvtab)[3L] <- type.measure

  out <- list(call = this.call, J = J, cv = cvlist, cvtab = cvtab,
              J.min = J[best.i], lambda.min = best$lambda.min,
              lambda.1se = best$lambda.1se, cvm.min = best$cvm.min,
              cvsd.min = best$cvsd.min, nzero.min = best$nzero.min,
              type.measure = type.measure, measure.name = best.name,
              nfolds = nfolds, foldid = foldid, wall.fit = best.spec)
  class(out) <- "cv.wall"

  return(out)

}

#' @rdname cv.wall
#' @param digits Significant digits in the printout.
#' @export
print.cv.wall <- function(x, digits = max(3L, getOption("digits") - 3L), ...){
  cat("\nCross-validated wall classifier\n\n")
  cat("Call:", deparse(x$call), "\n\n", sep = " ")
  cat("Measure:", x$measure.name, "(", x$nfolds, "folds )\n\n")
  tab <- x$cvtab
  tab$lambda.min <- signif(tab$lambda.min, digits)
  tab[[x$type.measure]] <- signif(tab[[x$type.measure]], digits)
  tab$sd <- signif(tab$sd, digits)
  print(tab, row.names = FALSE)
  cat("\nSelected: J =", x$J.min,
      "with lambda.min =", signif(x$lambda.min, digits),
      "( lambda.1se =", signif(x$lambda.1se, digits), ")\n")
  invisible(x)
}

#' @title Predictions and coefficients from a cross-validated wall classifier
#'
#' @description Obtains predictions and coefficients from a
#' \command{\link{cv.wall}} object, using the resolution level selected by
#' the cross-validation and, by default, the value \code{lambda.min} of the
#' penalty parameter.
#'
#' @param object A fitted object of class \code{"cv.wall"}.
#' @param newx Matrix (or data frame) of covariates at which predictions are
#'   to be made. See \command{\link{predict.wall}}.
#' @param s Value(s) of the penalty parameter at which predictions are
#'   computed: \code{"lambda.min"} (default, the minimizer of the
#'   cross-validated measure), \code{"lambda.1se"} or a numeric value.
#' @param ... Further arguments passed to \command{\link{predict.wall}}
#'   (e.g., \code{type}).
#'
#' @return See \command{\link{predict.wall}}.
#'
#' @seealso \command{\link{cv.wall}}, \command{\link{predict.wall}}
#'
#' @author Michel H. Montoril \email{michel@@ufscar.br}
#'
#' @examples
#' set.seed(1)
#' n <- 300
#' x <- matrix(runif(2*n), n, 2)
#' y <- rbinom(n, 1, 1/(1 + exp(-5*(x[, 1] - 0.5))))
#' cvfit <- cv.wall(x, y, J = 2:3, filter.size = 8, nfolds = 5)
#'
#' predict(cvfit, x[1:5, ], type = "response")
#' predict(cvfit, x[1:5, ], s = "lambda.1se", type = "class")
#' head(coef(cvfit))
#'
#' @keywords classif
#' @method predict cv.wall
#' @export
predict.cv.wall <- function(object, newx, s = c("lambda.min", "lambda.1se"),
                            ...){
  if(is.character(s)){
    s <- match.arg(s)
    s <- object[[s]]
  }
  predict(object$wall.fit, newx = newx, s = s, ...)
}

#' @rdname predict.cv.wall
#' @method coef cv.wall
#' @export
coef.cv.wall <- function(object, s = c("lambda.min", "lambda.1se"), ...){
  if(is.character(s)){
    s <- match.arg(s)
    s <- object[[s]]
  }
  coef(object$wall.fit, s = s, ...)
}

#' @title Plot the cross-validation curves of a cv.wall object
#'
#' @description Plots the cross-validated measure as a function of
#' \eqn{\log(\lambda)}, with one curve per candidate resolution level
#' \code{J}. The error bars (one standard error) are drawn for the selected
#' \code{J}, and the selected values \code{lambda.min} and \code{lambda.1se}
#' are marked by vertical dotted lines.
#'
#' @param x A fitted object of class \code{"cv.wall"}.
#' @param se Logical. Should the one-standard-error bars of the selected
#'   \code{J} be drawn? Default is \code{TRUE}.
#' @param col Vector of colors, one per candidate \code{J}.
#' @param legend.pos Position of the legend, as in
#'   \command{\link[graphics]{legend}}, or \code{NULL} to omit it.
#' @param ... Further graphical parameters passed to
#'   \command{\link[graphics]{plot}}.
#'
#' @return Invisibly, the object \code{x}.
#'
#' @seealso \command{\link{cv.wall}}
#'
#' @author Michel H. Montoril \email{michel@@ufscar.br}
#'
#' @examples
#' set.seed(1)
#' n <- 300
#' x <- matrix(runif(2*n), n, 2)
#' y <- rbinom(n, 1, 1/(1 + exp(-5*(x[, 1] - 0.5))))
#' cvfit <- cv.wall(x, y, J = 1:3, filter.size = 8, nfolds = 5)
#' plot(cvfit)
#'
#' @keywords classif hplot
#' @method plot cv.wall
#' @importFrom graphics lines points segments legend abline
#' @export
plot.cv.wall <- function(x, se = TRUE, col = NULL,
                         legend.pos = "topleft", ...){

  nJ <- length(x$J)
  if(is.null(col))
    col <- seq_len(nJ) + 1L
  col <- rep_len(col, nJ)

  logl <- lapply(x$cv, function(z) log(z$lambda))
  best <- x$cv[[match(x$J.min, x$J)]]

  ylim <- range(unlist(lapply(x$cv, `[[`, "cvm")))
  if(se)
    ylim <- range(ylim, best$cvup, best$cvlo)

  plot(NA, xlim = range(unlist(logl)), ylim = ylim,
       xlab = expression(log(lambda)), ylab = x$measure.name, ...)

  if(se){
    i <- match(x$J.min, x$J)
    segments(logl[[i]], best$cvlo, logl[[i]], best$cvup,
             col = "grey", lwd = 0.5)
  }
  for(i in seq_len(nJ)){
    lines(logl[[i]], x$cv[[i]]$cvm, col = col[i], lwd = 1.5)
    points(log(x$cv[[i]]$lambda.min), x$cv[[i]]$cvm.min, col = col[i],
           pch = 19)
  }
  abline(v = log(x$lambda.min), lty = 3)
  abline(v = log(x$lambda.1se), lty = 3)

  if(!is.null(legend.pos))
    legend(legend.pos, legend = paste("J =", x$J), col = col, lty = 1,
           lwd = 1.5, bty = "n")

  invisible(x)

}
