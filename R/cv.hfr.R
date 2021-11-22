#' @name cv.hfr
#' @title Cross validation for a Hierarchical Feature Regression
#' @description The hierarchical feature regression is fitted by estimating a semi-supervised
#' hierarchical graph based on bi-variate partial correlations, decomposing a least
#' squares estimate along the estimated hierarchy and shrinking coefficients along the branches of the tree.
#'
#' @details This function fits an HFR to a grid of penalty or factor values. The results is a
#' matrix of coefficients with one column for each hyperparameter. By evaluating all hyperparameters
#' in a single function, the speed of the algorithm can be improved substantially (e.g. by estimating
#' level-specific regressions only once).
#'
#' @param x Input matrix, of dimension \eqn{(N\times p)}{(N x p)}; each row is an observation vector.
#' @param y Response variable.
#' @param penalty_grid A vector of penalties on the effective degrees of freedom of the regression.
#' @param factors_grid A vector of target effective degrees of freedom of the regression.
#' @param q The quantile cut-off (in terms of information contributed) above which to consider levels in the hierarchy.
#' @param intercept Should intercept be fitted (default=TRUE).
#' @param standardize Logical flag for x variable standardization prior to fitting the model. The coefficients are always returned on the original scale. Default is \code{standardize=TRUE}.
#' @param cluster_method hierarchical cluster algorithm used to construct an asset hierarchy.
#' @return An 'hfr' regression object.
#' @author Johann Pfitzinger
#' @references
#'
#' @examples
#' x = matrix(rnorm(100 * 20), 100, 20)
#' y = rnorm(100)
#' fit = cv.hfr(x, y, factors_grid = seq(0, 1, by = 0.1))
#' coef(fit)
#'
#' @export
#'
#' @seealso \code{hfr}, \code{coef} and \code{predict} methods
#'
#' @importFrom quadprog solve.QP


cv.hfr <- function(
  x,
  y,
  penalty_grid = NULL,
  factors_grid = NULL,
  q = NULL,
  intercept = TRUE,
  standardize = TRUE,
  cluster_method = c("DIANA", "single", "complete", "average", "ward")
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

  if (is.null(penalty_grid) & is.null(factors_grid)) {
    warning("Both 'penalty' and 'factors' are zero. Setting 'penalty = 0'")
    penalty_grid <- 0
  }

  if (!is.null(factors_grid)) {
    if (any(factors_grid > 1) || any(factors_grid < 0)) {
      stop("'factors' must be between 0 and 1.")
    }
  }

  if (!is.null(penalty_grid)) {
    grid_size <- length(penalty_grid)
  } else {
    grid_size <- length(factors_grid)
  }

  # Get feature names
  var_names <- colnames(x)
  if (is.null(var_names)) var_names <- paste("X", 1:ncol(x), sep = ".")
  if (intercept) var_names <- c("intercept", var_names)

  if (is.null(q)) {
    q <- min(nvars, sqrt(nobs)) / nvars
  }

  if (any(apply(x, 2, sd)==0))
    stop("Features can not have a standard deviation of zero.")

  if (standardize) {
    standard_mean <- apply(x, 2, mean)
    standard_sd <- apply(x, 2, sd)
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

  beta_mat <- c()
  for (i in 1:grid_size) {

    if (!is.null(penalty_grid)) {
      penalty_term <- -v$dof * penalty_grid[i] * nobs
      opt <- solve.QP(Dmat = Dmat,
                      dvec = dvec + penalty_term,
                      Amat = Amat,
                      bvec = bvec)
    } else {
      dof_constraint <- 1 + factors_grid[i] * (nvars-1-1e-8)
      opt <- solve.QP(Dmat = Dmat,
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

    beta_mat <- cbind(beta_mat, beta)

  }

  if (!is.null(penalty_grid)) {
    colnames(beta_mat) <- penalty_grid
  } else {
    colnames(beta_mat) <- factors_grid
  }


  if (intercept) fitted <- cbind(1, x) %*% beta_mat else fitted <- x %*% beta_mat
  resid <- y - fitted

  out <- list(
    call = match.call(),
    coefficients = beta_mat,
    factors_grid = factors_grid,
    penalty_grid = penalty_grid,
    fitted.values = fitted,
    residuals = resid,
    x = x,
    y = y,
    df = v$dof,
    intercept = intercept,
    BIC = nobs * log(mean(resid^2)) + 2 * log(nobs) * sum(v$dof * opt.par)
  )

  class(out) <- "cv.hfr"

  return(out)

}
