#' @name hfr
#' @title Fit a hierarchical feature regression
#' @description HFR is a regularized regression estimator that decomposes a least squares
#' regression along a supervised hierarchical graph, and shrinks the edges of the
#' estimated graph to regularize parameters. The algorithm leads to group shrinkage in the
#' regression parameters and a reduction in the effective model degrees of freedom.
#'
#' @details Shrinkage can be imposed by targeting an explicit effective degrees of freedom.
#' Setting the argument \code{kappa} to a value between \code{0} and \code{1} controls
#' the effective degrees of freedom of the fitted object as a percentage of \eqn{p}{p}.
#' When \eqn{p > N}{p > N} \code{kappa} is a percentage of \eqn{(N - 2)}{(N - 2)}.
#' If no \code{kappa} is set, a linear regression with \code{kappa = 1} is
#' estimated.
#'
#' Hierarchical clustering is performed using \code{hclust}. The default is set to
#' ward.D2 clustering but can be overridden by passing a method argument to \code{...}.
#'
#' For high-dimensional problems, the hierarchy becomes very large. Setting \code{q} to a value below 1
#' reduces the number of levels used in the hierarchy. \code{q} represents a quantile-cutoff of the amount of
#' variation contributed by the levels. The default (\code{q = NULL}) considers all levels.
#'
#' @param x Input matrix or data.frame, of dimension \eqn{(N\times p)}{(N x p)}; each row is an observation vector.
#' @param y Response variable.
#' @param weights an optional vector of weights to be used in the fitting process. Should be NULL or a numeric vector. If non-NULL, weighted least squares is used for the level-specific regressions.
#' @param kappa The target effective degrees of freedom of the regression as a percentage of \eqn{p}{p}.
#' @param q Thinning parameter representing the quantile cut-off (in terms of contributed variance) above which to consider levels in the hierarchy. This can used to reduce the number of levels in high-dimensional problems. Default is no thinning.
#' @param intercept Should intercept be fitted. Default is \code{intercept=TRUE}.
#' @param standardize Logical flag for x variable standardization prior to fitting the model. The coefficients are always returned on the original scale. Default is \code{standardize=TRUE}.
#' @param partial_method Indicate whether to use pairwise partial correlations, or shrinkage partial correlations.
#' @param ridge_lambda Optional penalty for level-specific regressions (useful in high-dimensional case)
#' @param ...  Additional arguments passed to \code{hclust}.
#' @return An 'hfr' regression object.
#' @author Johann Pfitzinger
#' @references
#' Pfitzinger, J. (2022).
#' Cluster Regularization via a Hierarchical Feature Regression.
#' arXiv 2107.04831[statML]
#'
#' @examples
#' x = matrix(rnorm(100 * 20), 100, 20)
#' y = rnorm(100)
#' fit = hfr(x, y, kappa = 0.5)
#' coef(fit)
#'
#' @export
#'
#' @seealso \code{\link{cv.hfr}}, \code{\link{se.avg}}, \code{\link{coef}}, \code{\link{plot}} and \code{\link{predict}} methods
#'
#' @importFrom stats sd


hfr <- function(
  x,
  y,
  weights = NULL,
  kappa = 1,
  q = NULL,
  intercept = TRUE,
  standardize = TRUE,
  partial_method = c("pairwise", "shrinkage"),
  ridge_lambda = 0,
  ...
  ) {

  if (is.null(nobs <- nrow(x)))
    stop("'x' must be a matrix")
  if (nobs == 0L)
    stop("0 (non-NA) cases")
  nvars <- ncol(x)
  if (nvars == 0L) {
    return(list(coefficients = numeric(), residuals = y,
                fitted.values = 0 * y, dof = 0, clust = NULL,
                intercept = intercept))
  }
  ny <- NCOL(y)
  if (is.matrix(y) && ny == 1)
    y <- drop(y)
  if (ny > 1)
    stop("'y' must be a single response variable")
  if (NROW(y) != nobs)
    stop("incompatible dimensions")

  if (any(is.na(y)) || any(is.na(x)))
    stop("'NA' values in 'x' or 'y'")

  if (kappa > 1 || kappa < 0) {
    stop("'kappa' must be between 0 and 1")
  }

  if (!is.null(weights)) {
    if (length(weights) != nobs)
      stop("'weights' must have same length as 'y'")
    if (any(is.na(weights)))
      stop("'NA' values in 'weights'")
    if (any(weights < 0))
      stop("'weights' can only contain positive numerical values")
    wts <- sqrt(weights)
  } else {
    wts <- rep(1, nobs)
  }

  partial_method = match.arg(partial_method)

  # Get feature names
  var_names <- colnames(x)
  if (is.null(var_names)) var_names <- paste("X", 1:ncol(x), sep = ".")
  if (intercept) var_names <- c("(Intercept)", var_names)

  # Convert 'x' to matrix
  x <- data.matrix(x)

  if (is.null(q)) {
    q <- 1
  }

  if (any(apply(x, 2, stats::sd)==0))
    stop("Features can not have a standard deviation of zero.")

  if (standardize) {
    standard_mean <- apply(x, 2, mean)
    standard_sd <- apply(x, 2, stats::sd)
    if (intercept) {
      xs <- as.matrix(scale(x, center = standard_mean, scale = standard_sd))
    } else {
      xs <- as.matrix(scale(x, scale = standard_sd, center = FALSE))
    }
  } else {
    xs <- x
  }

  v = .get_level_reg(xs, y, wts, nvars, nobs, q, intercept, partial_method, ridge_lambda, ...)
  meta_opt <- .get_meta_opt(y, kappa, nvars, nobs, var_names, standardize, intercept, standard_sd, standard_mean, v)

  beta <- drop(meta_opt$beta)
  opt_par <- drop(meta_opt$opt_par)

  if (intercept) fitted <- as.numeric(cbind(1, x) %*% beta) else fitted <- as.numeric(x %*% beta)
  resid <- as.numeric(y - fitted)

  nlevels <- dim(v$fit_mat)[2]
  fit_mat_adj <- v$fit_mat
  theta <- cumsum(opt_par)
  for (i in 1:(nlevels-1)) fit_mat_adj[,i] <- fit_mat_adj[,i] - fit_mat_adj[,i+1]
  explained_variance <- sapply(nlevels:1, function(i) {
    if (i == nlevels) {
      fit <- fit_mat_adj[,i:nlevels] * theta[i:nlevels]
    } else {
      fit <- fit_mat_adj[,i:nlevels] %*% theta[i:nlevels]
    }
    return(1 - sum((fit - y)^2) / sum((y - mean(y))^2))
  })

  out <- list(
    call = match.call(),
    coefficients = beta,
    kappa = kappa,
    fitted.values = fitted,
    residuals = resid,
    x = x,
    y = y,
    df = round(as.numeric(crossprod(v$dof, opt_par)), 4) + intercept,
    hgraph = list(cluster_object = v$clust, shrinkage_vector = opt_par,
                  included_levels = v$included_levels,
                  explained_variance = explained_variance,
                  full_level_output = v),
    intercept = intercept
  )

  class(out) <- "hfr"

  return(out)

}
