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
#' \doi{10.1006/acha.1993.1005}.
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
#'   \command{\link{plot.wall}}, \command{\link{wbasis}},
#'   \command{\link{wtable}}, \command{\link[glmnet]{glmnet}}
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
#' @importFrom Matrix Matrix sparseMatrix drop0
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

#' @title Plot the regularization paths, the fitted components or the
#'   network structure of a wall classifier
#'
#' @description Three displays of a fitted \command{\link{wall}} object are
#' available. With \code{type = "path"}, the norm of the estimated
#' coefficients of each covariate is plotted against \eqn{\log(\lambda)},
#' summarizing the whole LASSO path and the order in which the covariates
#' enter the classifier. With \code{type = "components"}, the estimated
#' additive components \eqn{\hat{f}_l} are evaluated on a grid and plotted,
#' one panel per covariate, at a single value of the penalty parameter. With
#' \code{type = "network"}, the classifier is drawn as a layered graph, in
#' the fashion of a neural network diagram: the wavelet basis feeds the
#' additive components, which are summed into the log-odds, which the
#' logistic link turns into the predicted probability.
#'
#' @param x A fitted object of class \code{"wall"}.
#' @param type The display: \code{"path"} (default), the coefficient norms
#'   along the whole \eqn{\lambda} path; \code{"components"}, the fitted
#'   additive components at a given \eqn{\lambda}; or \code{"network"}, the
#'   layered structure of the classifier at a given \eqn{\lambda}.
#' @param s Value of the penalty parameter \eqn{\lambda}. Required by
#'   \code{type = "components"} and \code{type = "network"}, where it
#'   defaults to the smallest value of the path (the least penalized fit);
#'   typically, one passes here the value selected by
#'   \command{\link{cv.wall}}. With \code{type = "path"}, if \code{s} is
#'   provided, it is marked by a vertical dotted line and used to rank the
#'   covariates.
#' @param which Optional subset of covariates to be displayed, given as
#'   names, indices or a logical vector. By default (\code{which = NULL}) the
#'   covariates are selected automatically, as described below.
#' @param max.vars Maximum number of covariates displayed by the automatic
#'   selection: the ones with the largest coefficient norm. Default is
#'   \code{max.vars = 9}, which fills the default grid of panels of
#'   \code{type = "components"}. Use \code{max.vars = Inf} for all of them.
#'   Ignored if \code{which} is provided.
#' @param nonzero Logical, used by \code{type = "components"} and
#'   \code{type = "network"}. If \code{TRUE} (default), the automatic
#'   selection only considers the covariates whose component is not
#'   identically zero at \code{s}; in \code{type = "network"}, it also
#'   restricts the display to the resolution levels that survived the LASSO.
#' @param norm The norm used to summarize the block of coefficients of each
#'   covariate in \code{type = "path"}: \code{"l2"} (default) or \code{"l1"}.
#' @param center Logical, used by \code{type = "components"}. If \code{TRUE}
#'   (default), each component is centered at zero over the plotting grid.
#' @param n.grid Number of grid points at which the components are
#'   evaluated.
#' @param newx Optional single observation (a vector, or a matrix or data
#'   frame with one row and the covariates used in the fit), used by
#'   \code{type = "network"}. If provided, the nodes of the graph display the
#'   forward pass of that observation through the classifier, as described
#'   below.
#' @param col Vector of colors. For \code{type = "path"}, one per displayed
#'   covariate, defaulting to \code{2:(k+1)}, or to a qualitative palette
#'   when more than seven covariates are displayed; for
#'   \code{type = "components"}, one per displayed covariate, defaulting to
#'   black; for \code{type = "network"}, one per resolution level, from the
#'   scaling block (when present) to \eqn{W_{J-1}}, defaulting to grey for
#'   the scaling block and to a qualitative palette for the detail levels.
#' @param legend.pos Position of the legend of \code{type = "path"} and
#'   \code{type = "network"}, as in \command{\link[graphics]{legend}}, or
#'   \code{NULL} to omit it. Defaults to \code{"bottom"} (horizontally, over
#'   the free strip below the graph) with \code{type = "network"}.
#' @param mfrow Layout of the panels of \code{type = "components"}, as the
#'   homonymous parameter of \command{\link[graphics]{par}}. By default, the
#'   components are laid out in three columns, i.e., in a 3 by 3 grid with
#'   the default \code{max.vars}, with fewer rows when fewer components are
#'   displayed.
#' @param xlab,ylab,main,ylim Usual graphical parameters. Their defaults
#'   depend on \code{type}; with \code{type = "components"}, the axis labels
#'   are taken from the name of each covariate and \code{ylim} is shared by
#'   all the panels, so that the components are comparable, and with
#'   \code{type = "network"} the axes are omitted.
#' @param ... Further graphical parameters passed to
#'   \command{\link[graphics]{plot}}.
#'
#' @details
#' The coefficients of each covariate form a block of the whole coefficient
#' vector (the basis functions of its own expansion), so the path of a
#' covariate is summarized by the norm of its block, one curve per covariate.
#' The number of \emph{active} covariates, i.e., those with at least one
#' nonzero coefficient, is reported on the top axis.
#'
#' The components are evaluated on an equally spaced grid of the range of
#' each covariate observed in the fit, and are plotted on the original scale
#' of the covariate. They are expressed in the log-odds scale, and only their
#' sum is identifiable (see \command{\link{wall}}); with
#' \code{center = TRUE}, each of them is centered over the grid, which makes
#' the panels comparable and puts the (arbitrary) constants into the
#' intercept.
#'
#' Both displays remain readable when the number of covariates is large. In
#' \code{type = "path"}, only the \code{max.vars} most important covariates
#' are colored and named in the legend, while the remaining ones are drawn in
#' grey, so that the overall picture is kept. Importance is measured by the
#' coefficient norm at \code{s}, when a value is given, and by the average
#' norm along the path otherwise, which favors the covariates entering the
#' model early over those that only grow at the unpenalized end of the path.
#' In \code{type = "components"}, the covariates that are not displayed are
#' simply omitted from the layout of panels, with a message reporting how
#' many they are: the selection keeps the \code{max.vars} largest components
#' at \code{s}, measured by \eqn{\|\hat{f}_l\|_2 = \|\hat{d}_l\|_2} (the
#' basis is orthonormal), which are the ones that matter for the fitted
#' log-odds. By default, they fill a 3 by 3 grid of panels, in decreasing
#' order of importance; \code{max.vars} and \code{mfrow} change the number
#' of components displayed and their layout. In both cases, \code{which}
#' overrides the automatic selection.
#'
#' The display of \code{type = "network"} arranges the classifier in four
#' layers, read from right to left: the predicted probability
#' \eqn{\hat{\eta}(x)}; the log-odds \eqn{\hat{h}(x)}; the additive
#' components \eqn{\hat{f}_l}, one node per covariate; and the wavelet basis,
#' grouped by resolution level. A node of the leftmost layer collects the
#' coefficients of one resolution level of one covariate, i.e., the scaling
#' block \eqn{V_{j_0}} (when it is not dropped) or a detail block \eqn{W_j},
#' \eqn{j_0 \le j \le J_l - 1}; its color identifies the level, its size
#' grows with the number of coefficients that survived the LASSO there, and
#' the width of its edge is the norm \eqn{\|\hat{d}_{l,j}\|_2} of the block,
#' i.e., how much that level contributes to the component. With the default
#' \code{nonzero = TRUE}, only the levels with at least one nonzero
#' coefficient are drawn, so the picture displays the structure selected by
#' the LASSO. The remaining two layers carry no free parameters: the
#' components are summed into the log-odds, together with the unpenalized
#' intercept \eqn{\hat{\beta}_0} (drawn as a separate bias node), and the
#' logistic link maps the sum into the probability of the positive class.
#' The nodes are annotated in the notation of \command{\link{wall}}, with a
#' covariate identified by its index when the data carry no column names,
#' and by its name otherwise.
#'
#' With \code{newx}, the same graph displays the forward pass of a single
#' observation: each node of the basis layer is labeled with the detail of
#' the component at that level, \eqn{\sum_k \hat{d}_{jk}^l \psi_{jk}(x_l)},
#' each component node with \eqn{\hat{f}_l(x_l)}, and the last two nodes with
#' \eqn{\hat{h}(x)} and \eqn{\hat{\eta}(x)}; the edges leaving a component
#' are dashed when its contribution to the log-odds is negative. These values
#' are the ones actually used by \command{\link{predict.wall}} -- they are
#' not centered, so that they add up to \eqn{\hat{h}(x)} -- and the display
#' shows which covariates and which resolution levels drove the
#' classification of that observation.
#'
#' @return Invisibly, a list with the plotted quantities: for
#'   \code{type = "path"}, the components \code{lambda} and \code{norms} (a
#'   matrix with one row per covariate and one column per \eqn{\lambda});
#'   for \code{type = "components"}, the components \code{x} and
#'   \code{components} (matrices with one column per displayed covariate,
#'   holding the grid and the fitted values); for \code{type = "network"},
#'   the component \code{nodes} (a data frame with one row per node of the
#'   basis layer, holding its covariate, resolution level, label, number of
#'   nonzero coefficients, block norm and, with \code{newx}, its value),
#'   together with \code{f} (the components), \code{intercept}, \code{h} and
#'   \code{prob}, the last three being \code{NA} without \code{newx}. In
#'   every case, \code{which} gives the indices of the displayed covariates.
#'
#' @seealso \command{\link{wall}}, \command{\link{cv.wall}},
#'   \command{\link{predict.wall}}
#'
#' @author Michel H. Montoril \email{michel@@ufscar.br}
#'
#' @examples
#' set.seed(123)
#' n <- 400
#' x <- matrix(runif(6*n), n, 6)
#' h <- 3*sin(2*pi*x[, 1]) + 8*(x[, 2] - 0.5)   # only two active covariates
#' y <- rbinom(n, 1, 1/(1 + exp(-h)))
#'
#' fit <- wall(x, y, J = 3, filter.size = 8)
#' plot(fit)                                    # coefficient paths
#' plot(fit, type = "components", s = 0.02)     # fitted components
#' plot(fit, type = "network", s = 0.02)        # structure of the classifier
#'
#' # The forward pass of a single observation through the network
#' plot(fit, type = "network", s = 0.02, newx = x[1, ])
#'
#' # With the penalty parameter chosen by cross-validation
#' cvfit <- cv.wall(x, y, J = 2:3, filter.size = 8, nfolds = 5)
#' plot(cvfit$wall.fit, s = cvfit$lambda.min)
#' plot(cvfit$wall.fit, type = "components", s = cvfit$lambda.min)
#' plot(cvfit$wall.fit, type = "network", s = cvfit$lambda.min)
#'
#' @keywords classif hplot
#' @method plot wall
#' @importFrom graphics lines abline legend axis mtext par points segments
#' @importFrom graphics text strheight
#' @importFrom grDevices hcl.colors
#' @importFrom stats approx
#' @export
plot.wall <- function(x, type = c("path", "components", "network"), s = NULL,
                      which = NULL, max.vars = 9, nonzero = TRUE,
                      norm = c("l2", "l1"), center = TRUE, n.grid = 256,
                      newx = NULL, col = NULL, legend.pos = "topleft",
                      mfrow = NULL, xlab = NULL, ylab = NULL, main = NULL,
                      ylim = NULL, ...){

  type <- match.arg(type)
  norm <- match.arg(norm)

  if(!is.null(s) && (length(s) != 1L || !is.finite(s) || s < 0))
    stop("'s' must be a single non-negative value of the penalty parameter.")

  # The layered graph leaves its bottom strip free, where a horizontal
  # legend of the resolution levels fits without covering any node.
  if(type == "network" && missing(legend.pos))
    legend.pos <- "bottom"

  out <- switch(type,
                path =
                  .wall_plot_path(x, s = s, which = which,
                                  max.vars = max.vars, norm = norm, col = col,
                                  legend.pos = legend.pos, xlab = xlab,
                                  ylab = ylab, main = main, ylim = ylim, ...),
                components =
                  .wall_plot_components(x, s = s, which = which,
                                        max.vars = max.vars, nonzero = nonzero,
                                        center = center, n.grid = n.grid,
                                        col = col, mfrow = mfrow, xlab = xlab,
                                        ylab = ylab, main = main, ylim = ylim,
                                        ...),
                network =
                  .wall_plot_network(x, s = s, which = which,
                                     max.vars = max.vars, nonzero = nonzero,
                                     newx = newx, col = col,
                                     legend.pos = legend.pos, xlab = xlab,
                                     ylab = ylab, main = main, ylim = ylim,
                                     ...))

  invisible(c(list(type = type), out))

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
# with the same convention adopted in wdensity(). The column ranges can be
# supplied through 'ranges' (a 2 x d matrix, as returned by apply(x, 2,
# range)) to avoid their recomputation over repeated calls on the same data
# -- e.g., across the grid of J in cv.wall(), where only eps changes.
.wall_rescale_pars <- function(x, eps, rescale, boundary, ranges = NULL){
  d <- ncol(x)
  if(!rescale){
    if(min(x) < 0 || max(x) > 1)
      warning("Some covariates lie outside [0,1] and rescale = FALSE. The wavelet basis is only valid in the unit interval; the estimates may be unfeasible.")
    return(list(location = rep_len(0, d), scale = rep_len(1, d)))
  }
  if(is.null(ranges))
    ranges <- apply(x, 2L, range)
  location <- scale <- numeric(d)
  for(l in seq_len(d)){
    rx <- ranges[, l]
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

# Evaluates the wavelet basis of one (already rescaled) covariate directly
# in sparse (CSC) format, without the dense n x 2^J intermediate. Periodic
# boundary only. The result is the one of Matrix(B, sparse = TRUE) applied to
# the dense evaluation: same sparsity pattern to the bit (explicit zeros
# dropped), and the same values up to the rounding of the arithmetic, which
# the two loops reassociate differently where the compiler contracts a product
# and a sum into a fused multiply-add. The constant phi_{00} column is already
# removed when obj$drop.phi is TRUE.
.wall_wbasis_sparse <- function(u, obj, J){
  fam <- which(tolower(substring(obj$family, 1, 1)) == c("d", "s", "c", "o"))
  if(length(fam) == 0)
    stop("Unknown family. The families available are 'Daublets', 'Symmlets' and 'Coiflets'. You can also use 'Own', if you provide your wavelet.filter.")
  if(is.null(obj$wavelet.filter)){
    if(fam == 4)
      stop("Provide your own filter or choose an available family.")
    wf <- 0
  }
  else{
    fam <- 4L
    wf <- obj$wavelet.filter
  }
  wtab <- if(is.null(obj$wavelet.table)) NULL
          else .match_wavelet_table(obj$wavelet.table, fam, obj$filter.size,
                                    obj$wavelet.filter)
  out <- .Call("_WaveBased_C_WavBasisSparse", as.double(u),
               as.integer(obj$j0), as.integer(J), as.integer(fam),
               as.integer(obj$filter.size), as.integer(obj$prec.wavelet),
               as.double(wf), wtab, as.logical(obj$drop.phi))
  drop0(sparseMatrix(i = out$i, p = out$p, x = out$x,
                     dims = c(length(u), 2L^J - as.integer(obj$drop.phi)),
                     index1 = FALSE))
}

# Builds the stacked design matrix of basis functions, one block per
# covariate. With clip = TRUE (predictions), the rescaled covariates are
# truncated to the interval used in the fit. With sparse storage and the
# periodic boundary, each block is built directly in CSC format by
# .wall_wbasis_sparse(), which yields exactly the same matrix as the dense
# evaluation followed by Matrix(B, sparse = TRUE).
.wall_design <- function(x, obj, clip = FALSE){
  d <- ncol(x)
  blocks <- vector("list", d)
  sparse.direct <- obj$sparse && obj$boundary == "periodic"
  for(l in seq_len(d)){
    u <- (x[, l] - obj$location[l])/obj$scale[l]
    if(clip)
      u <- pmin(pmax(u, obj$eps[l]), 1 - obj$eps[l])
    if(sparse.direct){
      B <- .wall_wbasis_sparse(u, obj, obj$J[l])
      colnames(B) <- .wall_colnames(obj$xnames[l], obj$j0, obj$J[l], obj$drop.phi)
      blocks[[l]] <- B
      next
    }
    B <- .wall_wbasis(u, obj, obj$J[l])
    if(obj$drop.phi)
      B <- B[, -1L, drop = FALSE]
    colnames(B) <- .wall_colnames(obj$xnames[l], obj$j0, obj$J[l], obj$drop.phi)
    blocks[[l]] <- if(obj$sparse) Matrix(B, sparse = TRUE) else B
  }
  do.call(cbind, blocks)
}

# Column indices of the design (equivalently, of the coefficient vector
# without the intercept) belonging to each covariate: the basis functions of
# its own expansion, in the order used by .wall_design().
.wall_blocks <- function(obj){
  nphi <- if(obj$drop.phi) 0L else 2L^obj$j0
  sizes <- as.integer(nphi + 2^obj$J - 2^obj$j0)
  ends <- cumsum(sizes)
  Map(seq.int, ends - sizes + 1L, ends)
}

# Norm of the block of coefficients of each covariate, at every value of
# lambda: a matrix with one row per covariate and one column per lambda. The
# blocks of the (sparse) coefficient matrix are densified one at a time,
# which keeps the memory bounded by a single covariate.
.wall_block_norms <- function(beta, blocks, norm){
  out <- matrix(0, length(blocks), ncol(beta))
  for(l in seq_along(blocks)){
    b <- as.matrix(beta[blocks[[l]], , drop = FALSE])
    out[l, ] <- if(norm == "l1") colSums(abs(b)) else sqrt(colSums(b^2))
  }
  out
}

# Default layout of the panels of plot.wall(type = "components"): three
# components per row, which fills a 3 x 3 grid with the default max.vars,
# and keeps the panels with the same shape when there are fewer of them.
.wall_mfrow <- function(k){
  nc <- min(3L, k)
  c(ceiling(k/nc), nc)
}

# Chooses the covariates to be displayed by plot.wall(). A user-supplied
# 'which' (names, indices or a logical vector) is honored as given;
# otherwise the covariates with the largest 'score' are selected, and the
# remaining ones are returned separately, as the background of the plot.
.wall_pick <- function(which, score, xnames, max.vars, drop.zero = FALSE){

  d <- length(xnames)

  if(!is.null(which)){
    if(is.character(which)){
      sel <- match(which, xnames)
      if(anyNA(sel))
        stop("Unknown covariate(s) in 'which': ",
             paste(which[is.na(sel)], collapse = ", "), ".")
    }
    else if(is.logical(which)){
      if(length(which) != d)
        stop("A logical 'which' must have one entry per covariate (", d, ").")
      sel <- seq_len(d)[which]
    }
    else{
      sel <- as.integer(which)
      if(anyNA(sel) || any(sel < 1L) || any(sel > d))
        stop("'which' must index the ", d, " covariate(s) of the fit.")
    }
    if(length(sel) == 0L)
      stop("'which' selected no covariate.")
    return(list(sel = sel, bg = setdiff(seq_len(d), sel)))
  }

  if(length(max.vars) != 1L || is.na(max.vars) || max.vars < 1)
    stop("'max.vars' must be a single positive value.")

  cand <- seq_len(d)
  if(drop.zero)
    cand <- cand[score[cand] > 0]
  if(length(cand) == 0L)
    stop("Every component is identically zero at this value of 's'. Try a smaller 's', or nonzero = FALSE.")

  cand <- cand[order(score[cand], decreasing = TRUE)]
  sel <- if(length(cand) > max.vars) cand[seq_len(max.vars)] else cand

  list(sel = sel, bg = setdiff(seq_len(d), sel))

}

# Coefficient paths of plot.wall(type = "path").
.wall_plot_path <- function(obj, s, which, max.vars, norm, col, legend.pos,
                            xlab, ylab, main, ylim, ...){

  blocks <- .wall_blocks(obj)
  nrm <- .wall_block_norms(obj$glmnet.fit$beta, blocks, norm)
  loglam <- log(obj$lambda)

  # Importance of a covariate: its coefficient norm at 's', when a value is
  # given, and its average norm along the path otherwise -- which favors the
  # covariates entering early over those that only grow at the unpenalized
  # end of the path.
  score <- if(is.null(s))
             rowMeans(nrm)
           else
             .wall_block_norms(coef(obj, s = s)[-1L, , drop = FALSE], blocks,
                               norm)[, 1L]

  pick <- .wall_pick(which, score, obj$xnames, max.vars)
  sel <- pick$sel
  k <- length(sel)

  if(is.null(col))
    col <- if(k > 7L) hcl.colors(k, "Dark 3") else seq_len(k) + 1L
  col <- rep_len(col, k)

  if(is.null(xlab))
    xlab <- expression(log(lambda))
  if(is.null(ylab))
    ylab <- paste0(if(norm == "l1") "L1" else "L2", " coefficient norm")
  if(is.null(ylim))
    ylim <- range(0, nrm[c(sel, pick$bg), , drop = FALSE])

  plot(NA, xlim = range(loglam), ylim = ylim, xlab = xlab, ylab = ylab,
       main = main, ...)

  # The number of covariates with at least one nonzero coefficient, which
  # plays here the role of the degrees of freedom of plot.glmnet(). A path
  # with a single lambda (a user-supplied value) has nothing to interpolate,
  # and is drawn as points.
  single <- length(loglam) < 2L

  if(!single){
    atx <- pretty(loglam)
    nact <- approx(loglam, colSums(nrm > 0), xout = atx, rule = 2,
                   method = "constant", f = 0)$y
    axis(3L, at = atx, labels = nact, tcl = NA)
    if(is.null(main))
      mtext("active covariates", side = 3L, line = 2, cex = par("cex.lab"))
  }

  for(l in pick$bg)
    lines(loglam, nrm[l, ], col = "grey85", type = if(single) "p" else "l")
  for(i in seq_len(k))
    lines(loglam, nrm[sel[i], ], col = col[i], lwd = 1.5,
          type = if(single) "p" else "l")

  if(!is.null(s))
    abline(v = log(s), lty = 3)

  if(!is.null(legend.pos)){
    leg <- obj$xnames[sel]
    lcol <- col
    if(length(pick$bg) > 0L){
      leg <- c(leg, paste(length(pick$bg), "other(s)"))
      lcol <- c(lcol, "grey85")
    }
    legend(legend.pos, legend = leg, col = lcol, lty = 1, lwd = 1.5,
           bty = "n")
  }

  rownames(nrm) <- obj$xnames
  list(lambda = obj$lambda, norm = norm, norms = nrm, which = sel)

}

# Fitted additive components of plot.wall(type = "components").
.wall_plot_components <- function(obj, s, which, max.vars, nonzero, center,
                                  n.grid, col, mfrow, xlab, ylab, main, ylim,
                                  ...){

  if(length(n.grid) != 1L || !is.finite(n.grid) || n.grid < 2)
    stop("'n.grid' must be a single integer larger than one.")
  if(!is.null(mfrow)){
    if(length(mfrow) != 2L || any(!is.finite(mfrow)) || any(mfrow < 1))
      stop("'mfrow' must be a vector with the number of rows and columns of the layout of panels.")
    mfrow <- as.integer(mfrow)
  }
  if(is.null(s))
    s <- min(obj$lambda)

  beta <- as.numeric(coef(obj, s = s)[, 1L])[-1L]   # without the intercept
  blocks <- .wall_blocks(obj)

  # The basis is orthonormal, so the L2 norm of a component is the norm of
  # its coefficients, and no evaluation is needed to rank the covariates.
  score <- vapply(blocks, function(i) sqrt(sum(beta[i]^2)), 0)
  pick <- .wall_pick(which, score, obj$xnames, max.vars, drop.zero = nonzero)
  sel <- pick$sel
  k <- length(sel)

  if(is.null(which) && k < length(obj$J))
    message(sprintf("Showing %d of %d covariates (%d with a nonzero component at s = %s). See 'which' and 'max.vars'.",
                    k, length(obj$J), sum(score > 0), format(signif(s, 3))))

  n.grid <- as.integer(n.grid)
  grid <- fit <- matrix(NA_real_, n.grid, k,
                        dimnames = list(NULL, obj$xnames[sel]))

  for(i in seq_len(k)){
    l <- sel[i]
    u <- seq(obj$eps[l], 1 - obj$eps[l], length.out = n.grid)
    B <- .wall_wbasis(u, obj, obj$J[l])
    if(obj$drop.phi)
      B <- B[, -1L, drop = FALSE]
    f <- as.numeric(B %*% beta[blocks[[l]]])
    grid[, i] <- obj$location[l] + obj$scale[l]*u
    fit[, i] <- if(center) f - mean(f) else f
  }

  if(is.null(col))
    col <- 1L
  col <- rep_len(col, k)

  if(is.null(ylim)){
    ylim <- range(fit)
    if(diff(ylim) == 0)
      ylim <- ylim + c(-0.5, 0.5)
  }

  if(is.null(mfrow) && k > 1L)
    mfrow <- .wall_mfrow(k)

  if(!is.null(mfrow)){
    op <- par(mfrow = mfrow, mar = c(3.6, 3.6, 1, 1) + 0.1,
              mgp = c(2.2, 0.8, 0))
    on.exit(par(op))
  }

  for(i in seq_len(k)){
    plot(grid[, i], fit[, i], type = "l", col = col[i], lwd = 1.5,
         ylim = ylim, main = main,
         xlab = if(is.null(xlab)) obj$xnames[sel[i]] else xlab,
         ylab = if(is.null(ylab)) paste0("f(", obj$xnames[sel[i]], ")")
                else ylab, ...)
    abline(h = 0, lty = 3, col = "grey60")
  }

  list(s = s, x = grid, components = fit, which = sel)

}

# Resolution-level blocks of the coefficients of one covariate: the scaling
# block V_{j0} (when it is not dropped) and the detail blocks W_j, for
# j0 <= j <= J_l - 1, in the order used by .wall_design(). The indices are
# relative to the block of the covariate, as given by .wall_blocks().
.wall_groups <- function(obj, l){
  j <- obj$j0:(obj$J[l] - 1L)
  sizes <- as.integer(c(if(!obj$drop.phi) 2^obj$j0, 2^j))
  ends <- cumsum(sizes)
  list(idx = Map(seq.int, ends - sizes + 1L, ends),
       level = c(if(!obj$drop.phi) obj$j0, j),
       label = c(if(!obj$drop.phi) paste0("V", obj$j0), paste0("W", j)))
}

# Numbers displayed on the nodes: enough digits to be read, without the
# common padding that format() applies to a vector.
.wall_fmt <- function(v) vapply(v, function(z) format(signif(z, 3)), "")

# Subscript identifying a covariate in the annotation of the network: its
# index, when the covariates carry the names X1, ..., Xd that .wall_x()
# gives to an unnamed design, and its own name otherwise.
.wall_sub <- function(obj, sel){
  if(identical(obj$xnames, paste0("X", seq_along(obj$J))))
    as.character(sel)
  else
    obj$xnames[sel]
}

# Mathematical annotation of the resolution levels, from their labels: the
# scaling block V_{j0} and the detail blocks W_j.
.wall_key_expr <- function(keys){
  lev <- as.integer(substring(keys, 2L))
  lapply(seq_along(keys), function(i)
    if(substr(keys[i], 1L, 1L) == "V") bquote(V[.(lev[i])])
    else bquote(W[.(lev[i])]))
}

# Label of a node, with the value of a forward pass underneath it when there
# is one. The labels are plotmath expressions, one per node, and the two
# lines are drawn separately: text() would left-align the lines of a
# multi-line string, and plotmath has no line break anyway. With
# adj[2] = 0.5 the pair is centered on the node; otherwise it hangs below it.
.wall_label <- function(x, y, label, value, adj, cex, yr){
  x <- rep_len(x, length(label))
  y <- rep_len(y, length(label))
  for(i in seq_along(label)){
    lab <- as.expression(label[[i]])
    if(is.null(value))
      text(x[i], y[i], lab, adj = adj, cex = cex, xpd = NA)
    else if(adj[2L] == 0.5){
      # Centered on the node: the lines are anchored back to back, so that
      # the pair reads as one block, apart from the labels of its neighbors.
      text(x[i], y[i] + 0.25*cex*yr, lab, adj = c(adj[1L], 0), cex = cex,
           xpd = NA)
      text(x[i], y[i] - 0.25*cex*yr, value[i], adj = c(adj[1L], 1), cex = cex,
           xpd = NA)
    }
    else{
      # Hanging below the node, clear of the label above it, whose height
      # depends on the subscripts of its annotation.
      text(x[i], y[i], lab, adj = adj, cex = cex, xpd = NA)
      text(x[i], y[i] - strheight(lab, cex = cex) - 0.5*cex*yr, value[i],
           adj = adj, cex = cex, xpd = NA)
    }
  }
  invisible(NULL)
}

# Layered graph of plot.wall(type = "network").
.wall_plot_network <- function(obj, s, which, max.vars, nonzero, newx, col,
                               legend.pos, xlab, ylab, main, ylim, ...){

  if(is.null(s))
    s <- min(obj$lambda)

  cf <- as.numeric(coef(obj, s = s)[, 1L])
  b0 <- cf[1L]
  beta <- cf[-1L]                                   # without the intercept
  blocks <- .wall_blocks(obj)
  d <- length(blocks)

  # With an observation, the graph displays its forward pass; a vector is
  # read as the covariates of a single observation.
  act <- !is.null(newx)
  if(act){
    newx <- .wall_x(newx)
    if(d > 1L && ncol(newx) == 1L && nrow(newx) == d)
      newx <- t(newx)
    if(ncol(newx) != d)
      stop("'newx' must have ", d, " column(s), as the data used in the fit.")
    if(nrow(newx) != 1L)
      stop("'newx' must contain a single observation, whose forward pass through the network is displayed.")
  }

  # The basis is orthonormal, so the covariates are ranked by the norm of
  # their coefficients, as in type = "components".
  score <- vapply(blocks, function(i) sqrt(sum(beta[i]^2)), 0)
  pick <- .wall_pick(which, score, obj$xnames, max.vars, drop.zero = nonzero)
  sel <- pick$sel
  k <- length(sel)

  if(is.null(which) && k < d)
    message(sprintf("Showing %d of %d covariates (%d with a nonzero component at s = %s). See 'which' and 'max.vars'.",
                    k, d, sum(score > 0), format(signif(s, 3))))

  # One node per resolution level of each displayed covariate: how many of
  # its coefficients survived the LASSO (the size of the node), how large
  # they are as a block (the width of its edge) and, in a forward pass, the
  # detail of the component at that level (the value of the node).
  nodes <- vector("list", k)
  fx <- rep(NA_real_, k)
  for(i in seq_len(k)){
    l <- sel[i]
    g <- .wall_groups(obj, l)
    bl <- beta[blocks[[l]]]
    if(act){
      u <- (newx[1L, l] - obj$location[l])/obj$scale[l]
      u <- min(max(u, obj$eps[l]), 1 - obj$eps[l])
      B <- as.numeric(.wall_wbasis(u, obj, obj$J[l]))
      if(obj$drop.phi)
        B <- B[-1L]
    }
    nnz <- vapply(g$idx, function(ii) sum(bl[ii] != 0), 0L)
    nrm <- vapply(g$idx, function(ii) sqrt(sum(bl[ii]^2)), 0)
    val <- if(act) vapply(g$idx, function(ii) sum(B[ii]*bl[ii]), 0)
           else rep(NA_real_, length(g$idx))
    if(act)
      fx[i] <- sum(val)
    keep <- if(nonzero) nnz > 0L else rep(TRUE, length(nnz))
    nodes[[i]] <- data.frame(var = rep(obj$xnames[l], sum(keep)),
                             level = g$level[keep], label = g$label[keep],
                             nnz = nnz[keep], norm = nrm[keep],
                             value = val[keep])
  }

  h <- if(act) as.numeric(predict(obj, newx, s = s)) else NA_real_
  prob <- if(act) 1/(1 + exp(-h)) else NA_real_

  # Vertical layout: the levels of a covariate are stacked in their own
  # block, consecutive blocks are separated by a blank slot, and each
  # component is aligned with the center of the block that builds it.
  ng <- vapply(nodes, nrow, 0L)
  slots <- sum(pmax(ng, 1L)) + (k - 1L)
  ylev <- vector("list", k)
  yf <- numeric(k)
  cur <- slots
  for(i in seq_len(k)){
    m <- max(ng[i], 1L)
    yi <- cur - seq_len(m) + 1
    ylev[[i]] <- if(ng[i] > 0L) yi else numeric(0)
    yf[i] <- mean(yi)
    cur <- cur - m - 1
  }
  yh <- (slots + 1)/2
  ybias <- yh + 1.2
  yhead <- max(slots, ybias) + 0.9

  # The palette spans every level of the fit, from the coarsest to the
  # finest, so that a color means the same across the covariates (and across
  # plots of the same object). The unpenalized scaling block, which is a
  # different kind of object, is kept neutral.
  keys <- c(if(!obj$drop.phi) paste0("V", obj$j0),
            paste0("W", obj$j0:(max(obj$J[sel]) - 1L)))
  ndet <- length(keys) - as.integer(!obj$drop.phi)
  if(is.null(col))
    col <- c(if(!obj$drop.phi) "grey45",
             hcl.colors(max(ndet, 2L), "Dark 3")[seq_len(ndet)])
  col <- rep_len(col, length(keys))

  nnz.max <- max(c(unlist(lapply(nodes, `[[`, "nnz")), 1L))
  nrm.max <- max(c(unlist(lapply(nodes, `[[`, "norm")), 0))
  fw <- if(act) abs(fx) else score[sel]
  fw.max <- max(c(fw, 0))
  lwd.of <- function(v, vmax) 0.7 + 3.3*(if(vmax > 0) v/vmax else 0)

  if(is.null(ylim))
    ylim <- c(if(identical(legend.pos, "bottom")) -1 else 0.2, yhead + 0.4)

  op <- par(mar = c(1, 1, if(is.null(main)) 1 else 2.4, 1) + 0.1)
  on.exit(par(op))

  plot(NA, xlim = c(0.42, 4.32), ylim = ylim, axes = FALSE, main = main,
       xlab = if(is.null(xlab)) "" else xlab,
       ylab = if(is.null(ylab)) "" else ylab, ...)

  # Radius, in user coordinates, of a symbol drawn with cex = 1: the labels
  # are placed just outside each node, whatever its size.
  usr <- par("usr")
  xr <- 0.5*par("cin")[1L]*diff(usr[1:2])/par("pin")[1L]
  yr <- 0.5*par("cin")[2L]*diff(usr[3:4])/par("pin")[2L]

  text(1:4, yhead, c("wavelet basis", "components", "logit", "probability"),
       font = 2, cex = 0.85)

  # Edges, drawn first so that the nodes sit on top of them. From the basis
  # to the components, they are the blocks of coefficients estimated by the
  # LASSO; the remaining layers carry no free parameter, and a dashed edge
  # marks a component pulling the log-odds down.
  for(i in seq_len(k)){
    nd <- nodes[[i]]
    if(nrow(nd) == 0L)
      next
    segments(1, ylev[[i]], 2, yf[i], col = col[match(nd$label, keys)],
             lwd = lwd.of(nd$norm, nrm.max))
  }
  segments(2, yf, 3, yh, col = "grey30", lwd = lwd.of(fw, fw.max),
           lty = if(act) ifelse(fx < 0, 2L, 1L) else 1L)
  segments(3, ybias, 3, yh, col = "grey60", lty = if(b0 < 0) 2L else 1L)
  segments(3, yh, 4, yh, col = "grey30", lwd = 2)
  text(3.6, yh + 0.55*yr, expression(1/(1 + e^-hat(h))), cex = 0.75,
       col = "grey30", adj = c(0.5, 0))

  # Nodes of the basis layer, one per resolution level of each covariate,
  # labeled by the block they hold and by how many of its coefficients
  # survived the LASSO.
  for(i in seq_len(k)){
    nd <- nodes[[i]]
    if(nrow(nd) == 0L)
      next
    cex.i <- 1.4 + 2.6*nd$nnz/nnz.max
    points(rep(1, nrow(nd)), ylev[[i]], pch = 21, col = "black",
           bg = col[match(nd$label, keys)], cex = cex.i)
    ke <- .wall_key_expr(nd$label)
    .wall_label(1 - cex.i*xr - 0.5*xr, ylev[[i]],
                lapply(seq_len(nrow(nd)),
                       function(z) bquote(.(ke[[z]])*
                                          .(sprintf(" (%d)", nd$nnz[z])))),
                if(act) .wall_fmt(nd$value), adj = c(1, 0.5), cex = 0.72,
                yr = yr)
  }

  # The remaining layers: components, log-odds (with the unpenalized
  # intercept as a bias node) and the predicted probability.
  points(rep(2, k), yf, pch = 21, col = "black", bg = "grey95", cex = 3.2)
  .wall_label(2, yf - 2.1*yr,
              lapply(.wall_sub(obj, sel),
                     function(nm) bquote(hat(f)[.(nm)](x[.(nm)]))),
              if(act) .wall_fmt(fx), adj = c(0.5, 1), cex = 0.75, yr = yr)

  points(3, ybias, pch = 22, col = "black", bg = "grey95", cex = 2.2)
  .wall_label(3 + 1.6*xr, ybias,
              list(bquote(hat(beta)[0] == .(signif(b0, 3)))), NULL,
              adj = c(0, 0.5), cex = 0.75, yr = yr)

  points(3, yh, pch = 21, col = "black", bg = "grey95", cex = 3.6)
  .wall_label(3, yh - 2.3*yr, list(quote(hat(h)(x))),
              if(act) .wall_fmt(h), adj = c(0.5, 1), cex = 0.75, yr = yr)

  points(4, yh, pch = 21, col = "black", bg = "grey95", cex = 4)
  .wall_label(4, yh - 2.5*yr,
              list(bquote(hat(P)(Y == .(obj$classnames[2L])))),
              if(act) .wall_fmt(prob), adj = c(0.5, 1), cex = 0.75, yr = yr)

  shown <- keys %in% unlist(lapply(nodes, `[[`, "label"))
  if(!is.null(legend.pos) && any(shown))
    legend(legend.pos, legend = as.expression(.wall_key_expr(keys[shown])),
           col = col[shown], pch = 19, pt.cex = 1.4, cex = 0.85, bty = "n",
           title = "resolution level",
           horiz = identical(legend.pos, "bottom"))

  nodes <- do.call(rbind, nodes)
  f <- if(act) fx else score[sel]
  names(f) <- obj$xnames[sel]
  list(s = s, nodes = nodes, f = f, intercept = b0, h = h, prob = prob,
       which = sel)

}
