#' @name hfr
#' @title Fit a Hierarchical Feature Regression
#' @description The hierarchical feature regression is fitted by estimating a semi-supervised
#' hierarchical graph based on bi-variate partial correlations, decomposing a least
#' squares estimate along the estimated hierarchy and shrinking coefficients along the branches of the tree.
#'
#' @details Shrinkage can be imposed using two mechanisms: using a penalty on the effective degrees of
#' freedom of the regression, or targetting an explicit effective degrees of freedom.
#' Setting the argument \code{penalty} to a positive value, implements the former approach, while
#' setting the argument \code{factors} to a value between 0 and 1 implements the latter.
#' If neither \code{penalty} or \code{factors} is set, a linear regression with \code{penalty = 0} is
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
#' @param penalty A penalty on the effective degrees of freedom of the regression.
#' @param factors The target effective degrees of freedom of the regression as a percentage of nvars.
#' @param q The quantile cut-off (in terms of information contributed) above which to consider levels in the hierarchy.
#' @param intercept Should intercept be fitted (default=TRUE).
#' @param standardize Logical flag for x variable standardization prior to fitting the model. The coefficients are always returned on the original scale. Default is \code{standardize=TRUE}.
#' @param cluster_method hierarchical cluster algorithm used to construct an asset hierarchy.
#' @return An 'hfr' regression object.
#' @author Johann Pfitzinger
#' @references
#' Pfitzinger, J. (2021).
#' Cluster Regularization via a Hierarchical Feature Regression.
#' _arXiv [statML] 2107.04831_
#'
#' @examples
#' x = matrix(rnorm(100 * 20), 100, 20)
#' y = rnorm(100)
#' fit = hfr(x, y, factors = 0.5)
#' coef(fit)
#'
#' @export
#'
#' @seealso \code{cv.hfr}, \code{coef}, \code{plot} and \code{predict} methods
#'
#' @importFrom quadprog solve.QP
#' @importFrom stats sd


hfr <- function(
  x,
  y,
  penalty = NULL,
  factors = NULL,
  q = NULL,
  intercept = TRUE,
  standardize = TRUE,
  cluster_method = c("complete", "single", "average", "ward.D2", "mcquitty", "median", "centroid")
  ) {

  cluster_method = match.arg(cluster_method)

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

  if (is.null(penalty) & is.null(factors)) {
    warning("both 'penalty' and 'factors' are zero, setting 'penalty = 0'")
    penalty <- 0
  }

  if (!is.null(factors)) {
    if (factors > 1 || factors < 0) {
      stop("'factors' must be between 0 and 1")
    }
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
  }

  v = .get_level_reg(xs, y, nvars, nobs, cluster_method, q, intercept)

  Dmat <- crossprod(v$fit_mat) * 2
  diag(Dmat) <- diag(Dmat) + 1e-8
  dvec <- 2*t(v$fit_mat) %*% y
  Amat <- diag(length(v$dof))
  Amat[upper.tri(Amat)] <- 1
  Amat <- cbind(Amat, -Amat)
  bvec <- c(rep(0, length(v$dof)), -rep(1, length(v$dof)))

  if (!is.null(penalty)) {
    penalty_term <- -v$dof * penalty * nobs
    opt <- quadprog::solve.QP(Dmat = Dmat,
                              dvec = dvec + penalty_term,
                              Amat = Amat,
                              bvec = bvec)
  } else {
    dof_constraint <- 1 + factors * (nvars-1-1e-8)
    opt <- quadprog::solve.QP(Dmat = Dmat,
                              dvec = dvec,
                              Amat = cbind(v$dof, Amat),
                              bvec = c(dof_constraint, bvec),
                              meq = 1)
  }

  opt.par <- opt$solution

  beta <- rowSums(t(t(v$coef_mat) * opt.par))
  names(beta) <- var_names

  # Rescale beta
  if (standardize) {
    if (intercept) {
      beta[-1] <- beta[-1] / standard_sd
      beta[1] <- beta[1] - crossprod(beta[-1], standard_mean)
    } else {
      beta <- beta / standard_sd
    }
  }

  if (intercept) fitted <- as.numeric(cbind(1, x) %*% beta) else fitted <- as.numeric(x %*% beta)
  resid <- as.numeric(y - fitted)

  out <- list(
    call = match.call(),
    coefficients = beta,
    penalty = penalty,
    factors = factors,
    fitted.values = fitted,
    residuals = resid,
    x = x,
    y = y,
    df = as.numeric(crossprod(v$dof, opt.par)),
    cluster_model = list(cluster_object = v$clust, shrinkage_vector = opt.par,
                         included_levels = v$included_levels),
    intercept = intercept,
    BIC = nobs * log(mean(resid^2)) + 2 * log(nobs) * sum(v$dof * opt.par)
  )

  class(out) <- "hfr"

  return(out)

}
