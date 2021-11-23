#' @name hfr
#' @title Fit a hierarchical feature regression
#' @description HFR is a regularized regression estimator that decomposes a least squares
#' regression along a semi-supervised hierarchical graph, and shrinks coefficients
#' along the branches of the tree. The algorithm leads to group shrinkage in the
#' regression parameters and a reduction in the effective model degrees of freedom.
#'
#' @details Shrinkage can be imposed by targeting an explicit effective degrees of freedom.
#' Setting the argument \code{nu} to a value between 0 and 1 controls the effective degrees of
#' freedom of the fitted object as a percentage of \eqn{p}{p}. When \eqn{p > N}{p > N}
#' 'nu' is a percentage of \eqn{(N - 2)}{(N - 2)}.
#' If no \code{nu} is set, a linear regression with \code{nu = 1} is
#' estimated.
#'
#' Hierarchical clustering is performed using \code{hclust}. Default is complete-linkage agglomerative nesting.
#'
#' For high-dimensional problems, the hierarchy becomes very large. Setting \code{q} to a value below 1
#' reduces the number of levels used in the hierarchy. \code{q} represents a quantile-cutoff of amount of
#' information contributed by the levels. The default (\code{q = 1}) considers all levels.
#'
#' @param x Input matrix, of dimension \eqn{(N\times p)}{(N x p)}; each row is an observation vector.
#' @param y Response variable.
#' @param nu The target effective degrees of freedom of the regression as a percentage of nvars.
#' @param q The quantile cut-off (in terms of information contributed) above which to consider levels in the hierarchy.
#' @param intercept Should intercept be fitted (default=TRUE).
#' @param standardize Logical flag for x variable standardization prior to fitting the model. The coefficients are always returned on the original scale. Default is \code{standardize=TRUE}.
#' @param ...  Additional arguments passed to \code{hclust}.
#' @return An 'hfr' regression object.
#' @author Johann Pfitzinger
#' @references
#' Pfitzinger, J. (2021).
#' Cluster Regularization via a Hierarchical Feature Regression.
#' arXiv 2107.04831[statML]
#'
#' @examples
#' x = matrix(rnorm(100 * 20), 100, 20)
#' y = rnorm(100)
#' fit = hfr(x, y, nu = 0.5)
#' coef(fit)
#'
#' @export
#'
#' @seealso \code{cv.hfr}, \code{coef}, \code{plot} and \code{predict} methods
#'
#' @importFrom stats sd


hfr <- function(
  x,
  y,
  nu = 1,
  q = NULL,
  intercept = TRUE,
  standardize = TRUE,
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

  if (nu > 1 || nu < 0) {
    stop("'nu' must be between 0 and 1")
  }

  # Get feature names
  var_names <- colnames(x)
  if (is.null(var_names)) var_names <- paste("X", 1:ncol(x), sep = ".")
  if (intercept) var_names <- c("intercept", var_names)

  if (is.null(q)) {
    q <- min(nvars, sqrt(nobs)) / nvars
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

  v = .get_level_reg(xs, y, nvars, nobs, q, intercept, ...)
  meta_opt <- .get_meta_opt(y, nu, nvars, nobs, var_names, standardize, intercept, standard_sd, standard_mean, v)

  beta <- drop(meta_opt$beta)
  opt_par <- drop(meta_opt$opt_par)

  if (intercept) fitted <- as.numeric(cbind(1, x) %*% beta) else fitted <- as.numeric(x %*% beta)
  resid <- as.numeric(y - fitted)

  out <- list(
    call = match.call(),
    coefficients = beta,
    nu = nu,
    fitted.values = fitted,
    residuals = resid,
    x = x,
    y = y,
    df = round(as.numeric(crossprod(v$dof, opt_par)), 4),
    cluster_model = list(cluster_object = v$clust, shrinkage_vector = opt_par,
                         included_levels = v$included_levels),
    intercept = intercept
  )

  class(out) <- "hfr"

  return(out)

}
