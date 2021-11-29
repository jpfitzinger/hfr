#' @name cv.hfr
#' @title Cross validation for a hierarchical feature regression
#' @description HFR is a regularized regression estimator that decomposes a least squares
#' regression along a semi-supervised hierarchical graph, and shrinks coefficients
#' along the branches of the tree. The algorithm leads to group shrinkage in the
#' regression parameters and a reduction in the effective model degrees of freedom.
#'
#' @details This function fits an HFR to a grid of 'nu' hyperparameter values. The result is a
#' matrix of coefficients with one column for each hyperparameter. By evaluating all hyperparameters
#' in a single function, the speed of the algorithm can be improved substantially (e.g. by estimating
#' level-specific regressions only once).
#'
#' When 'nfolds > 1', a cross validation is performed with shuffled data. Alternatively,
#' test slices can be passed to the function using the 'foldid' argument. The result
#' of the cross validation is given by 'best_nu' in the output object.
#'
#' @param x Input matrix, of dimension \eqn{(N\times p)}{(N x p)}; each row is an observation vector.
#' @param y Response variable.
#' @param nu_grid A vector of target effective degrees of freedom of the regression.
#' @param q The quantile cut-off (in terms of information contributed) above which to consider levels in the hierarchy.
#' @param intercept Should intercept be fitted (default=TRUE).
#' @param standardize Logical flag for x variable standardization prior to fitting the model. The coefficients are always returned on the original scale. Default is \code{standardize=TRUE}.
#' @param nfolds The number of folds for k-fold cross validation (default=10).
#' @param foldid n optional vector of values between 1 and nfolds identifying what fold each observation is in. If supplied, nfolds can be missing.
#' @param ...  Additional arguments passed to \code{hclust}.
#' @return A 'cv.hfr' regression object.
#' @author Johann Pfitzinger
#' @references
#' Pfitzinger, J. (2021).
#' Cluster Regularization via a Hierarchical Feature Regression.
#' arXiv 2107.04831[statML]
#'
#' @examples
#' x = matrix(rnorm(100 * 20), 100, 20)
#' y = rnorm(100)
#' fit = cv.hfr(x, y, nu_grid = seq(0, 1, by = 0.1))
#' coef(fit)
#'
#' @export
#'
#' @seealso \code{hfr}, \code{coef} and \code{predict} methods
#'
#' @importFrom quadprog solve.QP
#' @importFrom stats sd


cv.hfr <- function(
  x,
  y,
  nu_grid = seq(0, 1, by = 0.1),
  q = NULL,
  intercept = TRUE,
  standardize = TRUE,
  nfolds = 10,
  foldid = NULL,
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

  if (any(nu_grid > 1) || any(nu_grid < 0)) {
    stop("each 'nu' must be between 0 and 1.")
  }

  if (is.null(foldid))
    foldid = sample(rep(seq(nfolds), length = nobs))
  else nfolds = max(foldid)

  # Get feature names
  var_names <- colnames(x)
  if (is.null(var_names)) var_names <- paste("X", 1:ncol(x), sep = ".")
  if (intercept) var_names <- c("intercept", var_names)

  if (is.null(q)) {
    q <- min(nvars, sqrt(nobs)) / nvars
  }

  if (nfolds > 1) {
    mse <- c()
    for (i in 1:nfolds) {

      ix <- foldid == i
      x_pred <- x[ix,,drop=FALSE]
      y_pred <- y[ix]
      x_fit <- x[!ix,,drop=FALSE]
      y_fit <- y[!ix]

      if (any(apply(x_fit, 2, stats::sd)==0))
        stop("Features can not have a standard deviation of zero.")

      if (standardize) {
        standard_mean <- apply(x_fit, 2, mean)
        standard_sd <- apply(x_fit, 2, stats::sd)
        if (intercept) {
          xs <- as.matrix(scale(x_fit, center = standard_mean, scale = standard_sd))
        } else {
          xs <- as.matrix(scale(x_fit, scale = standard_sd, center = FALSE))
        }
      } else {
        xs = x_fit
      }

      v = .get_level_reg(xs, y_fit, nvars, nobs, q, intercept, ...)
      meta_opt <- .get_meta_opt(y_fit, nu_grid, nvars, nobs, var_names, standardize, intercept, standard_sd, standard_mean, v)

      beta_mat <- meta_opt$beta
      opt_par_mat <- meta_opt$opt_par

      if (intercept) pred <- cbind(1, x_pred) %*% beta_mat else pred <- x_pred %*% beta_mat
      mse <- rbind(mse, colMeans((y_pred - pred)^2))

    }
    cv_mse <- colMeans(mse)
    best_nu <- nu_grid[which.min(cv_mse)]
  } else {
    best_nu <- NULL
    cv_mse <- NULL
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
  meta_opt <- .get_meta_opt(y, nu_grid, nvars, nobs, var_names, standardize, intercept, standard_sd, standard_mean, v)

  beta_mat <- meta_opt$beta
  opt_par_mat <- meta_opt$opt_par

  if (intercept) fitted <- cbind(1, x) %*% beta_mat else fitted <- x %*% beta_mat
  resid <- y - fitted

  TSS <- sum(y^2)
  explained_variance <- 1 - apply(sweep(v$fit_mat, 1, y), 2, function(f) sum(f^2) / TSS)

  out <- list(
    call = match.call(),
    coefficients = beta_mat,
    nu_grid = nu_grid,
    best_nu = best_nu,
    cv_mse = cv_mse,
    fitted.values = fitted,
    residuals = resid,
    x = x,
    y = y,
    df = round(as.numeric(v$dof %*% opt_par_mat), 4) + intercept,
    hgraph = list(cluster_object = v$clust, shrinkage_vector = opt_par_mat,
                  included_levels = v$included_levels,
                  explained_variance = explained_variance,
                  full_level_output = v),
    intercept = intercept
  )

  class(out) <- "cv.hfr"

  return(out)

}
