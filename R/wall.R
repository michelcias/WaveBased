#' @title Wavelet-based Additive Logistic LASSO classifier
#'
#' @description The \code{wall} function fits a nonparametric binary
#' classifier based on an additive logistic regression model whose components
#' are expanded in orthonormal wavelet bases. The coefficients of the
#' expansion are estimated and regularized simultaneously by the LASSO,
#' through the \pkg{glmnet} package, for a fixed resolution level \code{J} and
#' a path of values of the penalty parameter \eqn{\lambda}. See
#' \command{\link{cv.wall}} for the cross-validated choice of \code{J} and
#' \eqn{\lambda}.
#'
#' @param x Matrix (or data frame) of covariates, with \eqn{n} rows
#'   (observations) and \eqn{d} columns (variables). A vector is interpreted
#'   as a single covariate. For the \code{print} method, \code{x} is an
#'   object of class \code{"wall"}.
#' @param y Response with two classes. It can be a factor with two levels, a
#'   logical vector or a numeric vector with two distinct values. The second
#'   level (in sorted order, for non-factors) is treated as the positive
#'   class, i.e., the class whose conditional probability is modeled.
#' @param J The resolution level of the multiresolution space onto which each
#'   additive component is projected. Either a single integer, used for all
#'   the covariates, or a vector with one level per covariate (advanced use).
#'   It must satisfy \code{J > j0}.
#' @param j0 The coarsest resolution level of the decomposed wavelet basis.
#'   Default is \code{j0 = 0}. With \code{boundary = "interval"}, larger
#'   values are required (an informative error reports the minimum).
#' @param family The family of wavelets to use, as in \command{\link{wbasis}}:
#'   \emph{Daublets}, \emph{Symmlets} or \emph{Coiflets} (case insensitive,
#'   first letter suffices), or \emph{Own} alongside \code{wavelet.filter}.
#' @param filter.size The size of the wavelet filter. The available values
#'   depend on the chosen family; see \command{\link{wbasis}}.
#' @param prec.wavelet The number of iterations of the Daubechies-Lagarias
#'   algorithm used to evaluate the basis functions at the data points.
#' @param wavelet.filter Use this to provide your own filter, together with
#'   \code{family = "Own"}. Do not use it, if you are not sure about what you
#'   are doing. See \command{\link{wbasis}}.
#' @param boundary The boundary treatment of the basis. Either
#'   \code{"periodic"} (default, the periodized basis) or \code{"interval"}
#'   (the boundary-corrected orthonormal basis of Cohen, Daubechies and Vial,
#'   1993). See \command{\link{wbasis}}.
#' @param rescale Logical. If \code{TRUE} (default), each covariate is
#'   linearly rescaled to the interval \eqn{[}\code{eps}\eqn{, 1-}\code{eps}%
#'   \eqn{]} before the basis evaluation, and the transformation is stored in
#'   the fitted object so that \command{\link{predict.wall}} applies it to new
#'   data. If \code{FALSE}, the covariates are assumed to lie in \eqn{[0,1]}.
#' @param eps A value in \eqn{[0, 0.5)} used by the rescaling above. By
#'   default, \code{eps} \eqn{= 1.9^{-J}} for the periodized basis, following
#'   Montoril et al. (2019), and \code{eps} \eqn{= 0} for
#'   \code{boundary = "interval"} (whose construction does not suffer from
#'   the wrap-around boundary effect). Only used if \code{rescale = TRUE}.
#' @param use.table Controls whether the basis functions are evaluated by
#'   table lookup with linear interpolation (see \command{\link{wtable}}),
#'   which is considerably faster than the exact Daubechies-Lagarias
#'   algorithm for large data sets. One of \code{"auto"} (default: a table is
#'   built when the estimated amount of evaluations pays off its construction
#'   cost), \code{"always"} or \code{"never"}. Not used with
#'   \code{boundary = "interval"} (whose evaluation is always exact).
#' @param wavelet.table An optional object created by \command{\link{wtable}},
#'   overriding \code{use.table}. Providing a table built with
#'   \code{cache = TRUE} avoids its reconstruction across calls and sessions.
#' @param sparse Controls whether the design matrix of basis functions is
#'   stored in the sparse format of the \pkg{Matrix} package before fitting,
#'   which speeds up \pkg{glmnet} considerably when the proportion of
#'   non-zero entries is small (large \code{J} and/or small
#'   \code{filter.size}). One of \code{"auto"} (default: sparse storage is
#'   used when the expected density of the design is small), \code{"always"}
#'   or \code{"never"}.
#' @param lambda Optional user-supplied sequence of values of the penalty
#'   parameter \eqn{\lambda}. By default, \pkg{glmnet} computes its own
#'   sequence; see \command{\link[glmnet]{glmnet}}.
#' @param standardize Logical. Should the columns of basis functions be
#'   internally standardized by \pkg{glmnet} before the fit? Default is
#'   \code{FALSE}, which corresponds exactly to the penalized estimator
#'   in Montoril (2026), since the wavelet basis is already orthonormal.
#' @param weights Optional vector of observation weights, passed to
#'   \command{\link[glmnet]{glmnet}}.
#' @param ... Further arguments passed to \command{\link[glmnet]{glmnet}}
#'   (e.g., \code{nlambda}, \code{lambda.min.ratio}, \code{alpha},
#'   \code{maxit}). For the \code{print} method, further arguments are
#'   ignored.
#'
#' @details
#' Let \eqn{(X, Y) \in [0,1]^d \times \{0, 1\}} and let
#' \eqn{\eta(x) = P(Y = 1 | X = x)}. The classifier is based on the additive
#' approximation of the log-odds (logit) function
#' \deqn{h(x) = \log\left(\frac{\eta(x)}{1 - \eta(x)}\right) \approx
#' \sum_{l=1}^{d} f_l(x_l),}
#' which overcomes the curse of dimensionality. Each component \eqn{f_l} is
#' approximated by its orthogonal projection onto the multiresolution space
#' \eqn{V_{J_l}}, whose (decomposed) basis is
#' \eqn{\{\phi_{j_0 k}\}_k \cup \{\psi_{jk}, j_0 \le j \le J_l - 1\}_k}, as
#' provided by \command{\link{wbasis}}. The coefficients of the resulting
#' additive logistic regression are estimated by the LASSO, i.e., by
#' minimizing the penalized negative log-likelihood
#' \deqn{-\frac{1}{n}\sum_{i=1}^{n}\left[y_i(\beta_0^* + z_i^\top d) -
#' \log(1 + e^{\beta_0^* + z_i^\top d})\right] + \lambda \|d\|_1,}
#' where \eqn{z_i} collects the wavelet basis functions evaluated at the
#' covariates of observation \eqn{i} and \eqn{d} collects their coefficients.
#' The \eqn{\ell_1} penalty acts only on the wavelet (detail) coefficients
#' \eqn{d_{jk}^l}, performing automatic resolution-level and variable
#' selection: the scaling function block of each covariate is left
#' unpenalized. In particular, with the default \code{j0 = 0} and the
#' periodized basis, the single scaling function \eqn{\phi_{00}} of each
#' covariate is constant equal to 1 on \eqn{[0,1]}; these columns are
#' dropped from the design and absorbed into the (unpenalized) intercept
#' \eqn{\beta_0^*}, which handles the identifiability of the sum of the
#' scaling coefficients, as discussed in Montoril (2026).
#'
#' Following Montoril (2026), when \code{rescale = TRUE} the covariates are
#' linearly mapped to \eqn{[\epsilon, 1 - \epsilon]}, with
#' \eqn{\epsilon = 1.9^{-J}} by default, keeping the data away from the
#' boundary strip where the periodization artifact concentrates. The fitted
#' object stores the transformation of each covariate, and
#' \command{\link{predict.wall}} reapplies it to new data, truncating the
#' rescaled values to \eqn{[\epsilon, 1 - \epsilon]} whenever new
#' observations fall outside the range of the training data.
#'
#' Individual components \eqn{f_l} are generally not identifiable without
#' further centering constraints (only their sum is); the classifier and its
#' predictions are not affected by this. See Montoril (2026) for details.
#'
#' For computational speed, the wavelet basis can be evaluated by table
#' lookup (\code{use.table}, \command{\link{wtable}}) and the design matrix
#' can be stored in sparse format (\code{sparse}). By default, both choices
#' are made automatically from the size of the problem. When a lookup table
#' is used (or provided), it is stored in the fitted object, so that
#' predictions on new data are also fast; set \code{use.table = "never"} for
#' lighter objects with exact basis evaluations.
#'
#' @return An object of class \code{"wall"}, i.e., a list with components
#'   \item{glmnet.fit}{The fitted \code{"glmnet"} object, containing the whole
#'   LASSO path; see \command{\link[glmnet]{glmnet}}.}
#'   \item{lambda}{The sequence of values of \eqn{\lambda} actually used.}
#'   \item{J, j0}{The resolution levels of the fit (\code{J} has one entry
#'   per covariate).}
#'   \item{boundary, family, filter.size, prec.wavelet, wavelet.filter}{The
#'   wavelet basis specification.}
#'   \item{wavelet.table}{The lookup table used to evaluate the basis, or
#'   \code{NULL} for exact evaluations.}
#'   \item{rescale, location, scale, eps}{The rescaling specification: each
#'   covariate is mapped as \eqn{(x_l - location_l)/scale_l}.}
#'   \item{drop.phi}{Logical, whether the constant scaling function columns
#'   were dropped (periodized basis with \code{j0 = 0}).}
#'   \item{sparse}{Logical, whether the design matrix was stored in sparse
#'   format.}
#'   \item{penalty.factor}{The penalty factors applied to the columns of the
#'   design matrix (0 for scaling functions, 1 for wavelets).}
#'   \item{nobs, nvars}{Number of observations and of basis functions.}
#'   \item{xnames, classnames}{Names of the covariates and labels of the two
#'   classes.}
#'   \item{call}{The matched call.}
#'
#' @references
#' Montoril, M. H. (2026). Wavelet-based classification: an approach via
#' additive logistic regression. \emph{Manuscript}.
#'
#' Cohen, A., Daubechies, I. and Vial, P. (1993). Wavelets on the Interval and
#' Fast Wavelet Transforms. \emph{Applied and Computational Harmonic
#' Analysis}, 1(1), 54--81,
#' \url{https://doi.org/10.1006/acha.1993.1005}.
#'
#' Friedman, J., Hastie, T. and Tibshirani, R. (2010). Regularization Paths
#' for Generalized Linear Models via Coordinate Descent. \emph{Journal of
#' Statistical Software}, 33(1), 1--22, \doi{10.18637/jss.v033.i01}.
#'
#' Montoril, M. H., Chang, W. and Vidakovic, B. (2019). Wavelet-Based
#' Estimation of Generalized Discriminant Functions. \emph{Sankhya B}, 81(2),
#' 318--349, \doi{10.1007/s13571-018-0158-1}.
#'
#' Tibshirani, R. (1996). Regression shrinkage and selection via the lasso.
#' \emph{Journal of the Royal Statistical Society: Series B}, 58(1),
#' 267--288.
#'
#' @seealso \command{\link{cv.wall}}, \command{\link{predict.wall}},
#'   \command{\link{wbasis}}, \command{\link{wtable}},
#'   \command{\link[glmnet]{glmnet}}
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
#' # Fit for a fixed resolution level (whole lambda path)
#' fit <- wall(x, y, J = 3, family = "Daublets", filter.size = 8)
#' fit
#'
#' # Predictions at a specific lambda
#' prob <- predict(fit, x, s = 0.01, type = "response")
#' yhat <- predict(fit, x, s = 0.01, type = "class")
#' mean(yhat != y)  # training misclassification rate
#'
#' @keywords classif
#' @importFrom glmnet glmnet
#' @importFrom Matrix Matrix
#' @export
wall <- function(x, y, J, j0 = 0, family = "Daublets", filter.size = 20,
                 prec.wavelet = 30, wavelet.filter,
                 boundary = c("periodic", "interval"),
                 rescale = TRUE, eps,
                 use.table = c("auto", "always", "never"), wavelet.table = NULL,
                 sparse = c("auto", "always", "never"),
                 lambda = NULL, standardize = FALSE, weights = NULL, ...){

  this.call <- match.call()
  boundary <- match.arg(tolower(boundary[1L]), c("periodic", "interval"))
  use.table <- match.arg(use.table)
  sparse <- match.arg(sparse)

  x <- .wall_x(x)
  n <- nrow(x)
  d <- ncol(x)

  resp <- .wall_response(y)
  if(length(resp$y) != n)
    stop("The number of observations in 'y' must match the number of rows of 'x'.")

  if(missing(J))
    stop("The resolution level 'J' must be provided. See cv.wall() for a cross-validated choice.")
  j0 <- .wall_j0(j0)
  J <- .wall_J(J, j0, d)

  wf.own <- if(missing(wavelet.filter)) NULL else wavelet.filter

  if(missing(eps))
    eps <- NULL
  eps <- .wall_eps(eps, J, j0, boundary, rescale)
  rs <- .wall_rescale_pars(x, eps, rescale, boundary)

  wtab <- .wall_table(use.table, wavelet.table, boundary, workload = n*d,
                      family = family, filter.size = filter.size,
                      prec.wavelet = prec.wavelet, wavelet.filter = wf.own)

  drop.phi <- boundary == "periodic" && j0 == 0L
  L <- if(is.null(wf.own)) filter.size else length(wf.own)

  obj <- list(J = J, j0 = j0, boundary = boundary, family = family,
              filter.size = filter.size, prec.wavelet = prec.wavelet,
              wavelet.filter = wf.own, wavelet.table = wtab,
              rescale = rescale, location = rs$location, scale = rs$scale,
              eps = eps, drop.phi = drop.phi,
              sparse = .wall_sparse(sparse, J, j0, L, drop.phi),
              xnames = colnames(x), classnames = resp$classnames, nobs = n)

  design <- .wall_design(x, obj, clip = FALSE)
  pf <- .wall_penalty(J, j0, drop.phi)

  fit <- glmnet::glmnet(x = design, y = resp$y, family = "binomial",
                        weights = weights, lambda = lambda,
                        penalty.factor = pf, standardize = standardize, ...)

  obj$glmnet.fit <- fit
  obj$lambda <- fit$lambda
  obj$penalty.factor <- pf
  obj$nvars <- ncol(design)
  obj$call <- this.call
  class(obj) <- "wall"

  return(obj)

}

#' @rdname wall
#' @param digits Significant digits in the printout.
#' @export
print.wall <- function(x, digits = max(3L, getOption("digits") - 3L), ...){
  cat("\nWavelet-based additive logistic LASSO classifier (wall)\n\n")
  cat("Call:", deparse(x$call), "\n\n", sep = " ")
  cat("  Wavelet basis  :",
      if(is.null(x$wavelet.filter))
        paste0(x$family, " (filter size ", x$filter.size, ")")
      else
        paste0("own filter (size ", length(x$wavelet.filter), ")"),
      "\n")
  cat("  Boundary       :", x$boundary, "\n")
  cat("  Levels         : j0 =", x$j0, "and J =",
      paste(unique(x$J), collapse = ", "), "\n")
  cat("  Classes        :", paste(x$classnames, collapse = " / "), "\n")
  cat("  Data           : n =", x$nobs, "observations, d =", length(x$J),
      "covariate(s)\n")
  cat("  Basis functions:", x$nvars,
      if(x$sparse) "(sparse design)" else "(dense design)", "\n")
  cat("  Lambda path    :", length(x$lambda), "values in [",
      signif(min(x$lambda), digits), ",", signif(max(x$lambda), digits),
      "]\n\n")
  invisible(x)
}

#' @title Predictions and coefficients from a wall classifier
#'
#' @description Obtains the estimated log-odds, class probabilities, class
#' labels or wavelet coefficients from a fitted \command{\link{wall}} object,
#' at the requested values of the penalty parameter \eqn{\lambda}.
#'
#' @param object A fitted object of class \code{"wall"}.
#' @param newx Matrix (or data frame) of covariates at which predictions are
#'   to be made, with the same columns used in the fit. A vector is
#'   interpreted as a single covariate.
#' @param s Value(s) of the penalty parameter \eqn{\lambda} at which
#'   predictions are computed. Default is \code{s = NULL}, which returns
#'   predictions for the whole sequence used in the fit. See
#'   \command{\link[glmnet]{predict.glmnet}}.
#' @param type Type of prediction: \code{"link"} returns the estimated
#'   log-odds \eqn{\hat{h}(x)}; \code{"response"} returns the estimated
#'   conditional probabilities \eqn{\hat{\eta}(x)} of the positive (second)
#'   class; and \code{"class"} returns the predicted class labels, i.e., the
#'   plug-in classifier \eqn{I(\hat{h}(x) > 0)} mapped back to the original
#'   labels of \code{y}.
#' @param ... Further arguments passed to
#'   \command{\link[glmnet]{predict.glmnet}} (e.g., \code{exact}).
#'
#' @details
#' The covariates in \code{newx} are rescaled with the transformation stored
#' in \code{object} (see \command{\link{wall}}) and truncated to the interval
#' \eqn{[\epsilon, 1-\epsilon]} used in the fit, so that new observations
#' beyond the range of the training data are safely handled. The wavelet
#' basis is evaluated with the same method used in the fit (lookup table or
#' exact algorithm).
#'
#' @return A matrix of predictions, with one row per observation of
#'   \code{newx} and one column per requested value of \eqn{\lambda}. For
#'   \code{type = "class"}, the matrix contains the class labels. For
#'   \code{coef.wall}, a sparse matrix of coefficients, whose row names
#'   identify the basis functions: e.g., \code{X2.psi3.5} is the wavelet
#'   \eqn{\psi_{3,5}} of the second covariate.
#'
#' @seealso \command{\link{wall}}, \command{\link{cv.wall}},
#'   \command{\link[glmnet]{predict.glmnet}}
#'
#' @author Michel H. Montoril \email{michel@@ufscar.br}
#'
#' @examples
#' set.seed(1)
#' n <- 300
#' x <- matrix(runif(2*n), n, 2)
#' y <- rbinom(n, 1, 1/(1 + exp(-5*(x[, 1] - 0.5))))
#' fit <- wall(x, y, J = 3, filter.size = 8)
#'
#' predict(fit, x[1:5, ], s = 0.05)                    # log-odds
#' predict(fit, x[1:5, ], s = 0.05, type = "response") # probabilities
#' predict(fit, x[1:5, ], s = 0.05, type = "class")    # labels
#' head(coef(fit, s = 0.05))
#'
#' @keywords classif
#' @method predict wall
#' @importFrom stats predict
#' @export
predict.wall <- function(object, newx, s = NULL,
                         type = c("link", "response", "class"), ...){

  type <- match.arg(type)

  if(missing(newx))
    stop("'newx' must be provided.")
  newx <- .wall_x(newx)
  if(ncol(newx) != length(object$J))
    stop("'newx' must have ", length(object$J), " column(s), as the data used in the fit.")

  design <- .wall_design(newx, object, clip = TRUE)

  pred <- predict(object$glmnet.fit, newx = design, s = s,
                  type = if(type == "class") "link" else type, ...)

  if(type == "class"){
    labels <- object$classnames[1L + (pred > 0)]
    dim(labels) <- dim(pred)
    dimnames(labels) <- dimnames(pred)
    return(labels)
  }

  return(pred)

}

#' @rdname predict.wall
#' @method coef wall
#' @importFrom stats coef
#' @export
coef.wall <- function(object, s = NULL, ...){
  predict(object$glmnet.fit, type = "coefficients", s = s, ...)
}


# ------------------------------------------------------------------------
# Internal helpers shared by wall() and cv.wall().
# ------------------------------------------------------------------------

# Coerces the covariates to a numeric matrix with column names.
.wall_x <- function(x){
  if(is.data.frame(x))
    x <- as.matrix(x)
  if(!is.matrix(x))
    x <- as.matrix(x)
  if(!is.numeric(x))
    stop("'x' must be numeric.")
  if(any(!is.finite(x)))
    stop("'x' must contain only finite values.")
  if(is.null(colnames(x)))
    colnames(x) <- paste0("X", seq_len(ncol(x)))
  x
}

# Coerces the two-class response to 0/1, keeping the original labels.
# The second class (second level, or largest value) is the positive one.
.wall_response <- function(y){
  if(any(is.na(y)))
    stop("'y' must not contain missing values.")
  if(is.factor(y)){
    if(nlevels(y) != 2L)
      stop("'y' must have exactly two classes.")
    classnames <- levels(y)
    y01 <- as.integer(y) - 1L
  }
  else{
    uy <- sort(unique(y))
    if(length(uy) != 2L)
      stop("'y' must have exactly two classes.")
    classnames <- as.character(uy)
    y01 <- as.integer(y == uy[2L])
  }
  list(y = y01, classnames = classnames)
}

# Validates the coarsest level.
.wall_j0 <- function(j0){
  if(length(j0) != 1L || !is.finite(j0) || j0 < 0 || j0 != round(j0))
    stop("'j0' must be a single non-negative integer.")
  as.integer(j0)
}

# Validates the resolution level(s) and recycles them to one per covariate.
.wall_J <- function(J, j0, d){
  if(!is.numeric(J) || any(!is.finite(J)) || any(J != round(J)))
    stop("'J' must contain only integer values.")
  if(length(J) != 1L && length(J) != d)
    stop("'J' must have length 1 or one entry per covariate (", d, ").")
  if(any(J <= j0))
    stop("'J' must be larger than 'j0' (= ", j0, ").")
  rep_len(as.integer(J), d)
}

# Resolves the boundary-avoidance constant of the rescaling, one entry per
# covariate: 1.9^(-J) for the periodized basis (Montoril et al., 2019) and 0
# for the interval basis or when the data are not rescaled.
.wall_eps <- function(eps, J, j0, boundary, rescale){
  if(!rescale || boundary == "interval")
    return(rep_len(0, length(J)))
  if(is.null(eps))
    return(1.9^(-J))
  if(!is.numeric(eps) || any(!is.finite(eps)) || any(eps < 0) || any(eps >= 0.5))
    stop("'eps' must belong to the interval [0, 0.5).")
  if(length(eps) != 1L && length(eps) != length(J))
    stop("'eps' must have length 1 or one entry per covariate (", length(J), ").")
  rep_len(eps, length(J))
}

# Location/scale of the linear map taking each covariate to [eps, 1-eps],
# with the same convention adopted in wdensity().
.wall_rescale_pars <- function(x, eps, rescale, boundary){
  d <- ncol(x)
  if(!rescale){
    if(min(x) < 0 || max(x) > 1)
      warning("Some covariates lie outside [0,1] and rescale = FALSE. The wavelet basis is only valid in the unit interval; the estimates may be unfeasible.")
    return(list(location = rep_len(0, d), scale = rep_len(1, d)))
  }
  location <- scale <- numeric(d)
  for(l in seq_len(d)){
    rx <- range(x[, l])
    if(rx[1L] == rx[2L])
      stop("Covariate ", colnames(x)[l], " is constant and cannot be rescaled.")
    a <- eps[l]*diff(rx)/(1 - 2*eps[l])
    location[l] <- rx[1L] - a
    scale[l] <- 2*a + diff(rx)
  }
  list(location = location, scale = scale)
}

# Resolves the lookup-table policy and builds the table if worthwhile. The
# "auto" rule compares the total amount of basis evaluations (workload, in
# points x covariates x designs) with the construction cost of the table,
# which grows with the filter size.
.wall_table <- function(use.table, wavelet.table, boundary, workload,
                        family, filter.size, prec.wavelet, wavelet.filter){
  if(boundary == "interval"){
    if(!is.null(wavelet.table))
      warning("wavelet.table is not used with boundary = \"interval\" and was ignored.")
    return(NULL)
  }
  if(!is.null(wavelet.table))
    return(wavelet.table)
  L <- if(is.null(wavelet.filter)) filter.size else length(wavelet.filter)
  build <- switch(use.table,
                  always = TRUE,
                  never = FALSE,
                  auto = workload >= 2000*L)
  if(!build)
    return(NULL)
  if(is.null(wavelet.filter))
    wtable(family = family, filter.size = filter.size,
           prec.wavelet = prec.wavelet, check = FALSE)
  else
    wtable(family = "Own", prec.wavelet = prec.wavelet,
           wavelet.filter = wavelet.filter, check = FALSE)
}

# Decides the sparse storage of the design. The expected density comes from
# the support of the basis: at level j, at most min(2^j, L-1) translates are
# non-zero at any given point.
.wall_sparse <- function(sparse, J, j0, L, drop.phi){
  if(sparse == "always")
    return(TRUE)
  if(sparse == "never")
    return(FALSE)
  p <- nnz <- 0
  for(Jl in J){
    p <- p + 2^Jl - if(drop.phi) 1L else 0L
    nnz <- nnz + (if(drop.phi) 0 else min(2^j0, L - 1)) +
      sum(pmin(2^(j0:(Jl - 1L)), L - 1))
  }
  p >= 128 && nnz/p < 0.4
}

# Penalty factors of the columns of the design: 0 for the scaling functions
# (unpenalized, they absorb the intercept term beta_0^*) and 1 for the
# wavelets, matching the l1 penalty on the detail coefficients only.
.wall_penalty <- function(J, j0, drop.phi){
  unlist(lapply(J, function(Jl){
    nphi <- if(drop.phi) 0L else 2^j0
    c(rep.int(0, nphi), rep.int(1, 2^Jl - 2^j0))
  }))
}

# Labels of the columns of one covariate block, e.g. "X1.psi2.3" for the
# wavelet psi_{2,3} of the first covariate.
.wall_colnames <- function(xname, j0, J, drop.phi){
  nm <- if(drop.phi) character(0L)
        else paste0(xname, ".phi", j0, ".", seq_len(2^j0) - 1L)
  for(j in j0:(J - 1L))
    nm <- c(nm, paste0(xname, ".psi", j, ".", seq_len(2^j) - 1L))
  nm
}

# Evaluates the wavelet basis of one (already rescaled) covariate, honoring
# the specification stored in the (possibly partial) wall object.
.wall_wbasis <- function(u, obj, J){
  if(is.null(obj$wavelet.filter))
    wbasis(u, j0 = obj$j0, J = J, family = obj$family,
           filter.size = obj$filter.size, prec.wavelet = obj$prec.wavelet,
           boundary = obj$boundary, wavelet.table = obj$wavelet.table)
  else
    wbasis(u, j0 = obj$j0, J = J, family = "Own",
           filter.size = obj$filter.size, prec.wavelet = obj$prec.wavelet,
           wavelet.filter = obj$wavelet.filter,
           boundary = obj$boundary, wavelet.table = obj$wavelet.table)
}

# Builds the stacked design matrix of basis functions, one block per
# covariate. With clip = TRUE (predictions), the rescaled covariates are
# truncated to the interval used in the fit.
.wall_design <- function(x, obj, clip = FALSE){
  d <- ncol(x)
  blocks <- vector("list", d)
  for(l in seq_len(d)){
    u <- (x[, l] - obj$location[l])/obj$scale[l]
    if(clip)
      u <- pmin(pmax(u, obj$eps[l]), 1 - obj$eps[l])
    B <- .wall_wbasis(u, obj, obj$J[l])
    if(obj$drop.phi)
      B <- B[, -1L, drop = FALSE]
    colnames(B) <- .wall_colnames(obj$xnames[l], obj$j0, obj$J[l], obj$drop.phi)
    blocks[[l]] <- if(obj$sparse) Matrix(B, sparse = TRUE) else B
  }
  do.call(cbind, blocks)
}
