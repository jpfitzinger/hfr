#' @name cv.hfr
#' @title Cross validation for a hierarchical feature regression
#' @description HFR is a regularized regression estimator that decomposes a least squares
#' regression along a supervised hierarchical graph, and shrinks the edges of the
#' estimated graph to regularize parameters. The algorithm leads to group shrinkage in the
#' regression parameters and a reduction in the effective model degrees of freedom.
#'
#' @details This function fits an HFR to a grid of \code{kappa} hyperparameter values. The result is a
#' matrix of coefficients with one column for each hyperparameter. By evaluating all hyperparameters
#' in a single function, the speed of the cross-validation procedure is improved substantially (since
#' level-specific regressions are estimated only once).
#'
#' When \code{nfolds > 1}, a cross validation is performed with shuffled data. Alternatively,
#' test slices can be passed to the function using the \code{foldid} argument. The result
#' of the cross validation is given by \code{best_kappa} in the output object.
#'
#' @param x Input matrix or data.frame, of dimension \eqn{(N\times p)}{(N x p)}; each row is an observation vector.
#' @param y Response variable.
#' @param weights an optional vector of weights to be used in the fitting process. Should be NULL or a numeric vector. If non-NULL, weighted least squares is used for the level-specific regressions.
#' @param kappa A vector of target effective degrees of freedom of the regression.
#' @param q Thinning parameter representing the quantile cut-off (in terms of contributed variance) above which to consider levels in the hierarchy. This can used to reduce the number of levels in high-dimensional problems. Default is no thinning.
#' @param intercept Should intercept be fitted. Default is \code{intercept=TRUE}.
#' @param standardize Logical flag for \code{x} variable standardization prior to fitting the model. The coefficients are always returned on the original scale. Default is \code{standardize=TRUE}.
#' @param nfolds The number of folds for k-fold cross validation. Default is \code{nfolds=10}.
#' @param foldid An optional vector of values between \code{1} and \code{nfolds} identifying what fold each observation is in. If supplied, \code{nfolds} can be missing.
#' @param partial_method Indicate whether to use pairwise partial correlations, or shrinkage partial correlations.
#' @param l2_penalty Optional penalty for level-specific regressions (useful in high-dimensional case)
#' @param ...  Additional arguments passed to \code{hclust}.
#' @return A 'cv.hfr' regression object.
#' @author Johann Pfitzinger
#' @references
#' Pfitzinger, Johann (2024). Cluster Regularization via a Hierarchical Feature Regression. _Econometrics and Statistics_ (in press). URL https://doi.org/10.1016/j.ecosta.2024.01.003.
#'
#' @examples
#' x = matrix(rnorm(100 * 20), 100, 20)
#' y = rnorm(100)
#' fit = cv.hfr(x, y, kappa = seq(0, 1, by = 0.1))
#' coef(fit)
#'
#' @export
#'
#' @seealso \code{\link{hfr}}, \code{\link{coef}}, \code{\link{plot}} and \code{\link{predict}} methods
#'
#' @importFrom quadprog solve.QP
#' @importFrom stats sd


cv.hfr <- function(
  x,
  y,
  weights = NULL,
  kappa = seq(0, 1, by = 0.1),
  q = NULL,
  intercept = TRUE,
  standardize = TRUE,
  nfolds = 10,
  foldid = NULL,
  partial_method = c("pairwise", "shrinkage"),
  l2_penalty = 0,
  ...
) {

  args = .check_args(x = x, y = y, intercept = intercept, kappa = kappa,
                     weights = weights, q = q, l2_penalty = l2_penalty, is_cv = TRUE)
  if (args$nvars == 0L) {
    return(list(coefficients = numeric(), residuals = y,
                fitted.values = 0 * y, dof = 0, clust = NULL,
                intercept = intercept))
  }
  partial_method = match.arg(partial_method)

  if (is.null(foldid))
    foldid = sample(rep(seq(nfolds), length = args$nobs))
  else nfolds = max(foldid)

  x <- args$x
  y <- args$y

  if (nfolds > 1) {
    mse <- c()
    for (i in 1:nfolds) {

      ix <- foldid == i
      x_pred <- x[ix,,drop=FALSE]
      y_pred <- y[ix]
      x_fit <- x[!ix,,drop=FALSE]
      y_fit <- y[!ix]

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

      v = .get_level_reg(xs, y_fit, args$wts[!ix], args$nvars, args$nobs, args$q,
                         intercept, partial_method, args$l2_penalty, ...)
      meta_opt <- .get_meta_opt(y_fit, args$kappa, args$nvars, args$nobs,
                                args$var_names_excl, standardize, intercept,
                                standard_sd, standard_mean, v)

      beta_mat <- meta_opt$beta
      opt_par_mat <- meta_opt$opt_par

      if (intercept) pred <- cbind(1, x_pred) %*% beta_mat else pred <- x_pred %*% beta_mat
      mse <- rbind(mse, colMeans((y_pred - pred)^2))

    }
    cv_mse <- colMeans(mse)
    best_kappa <- args$kappa[which.min(cv_mse)]
  } else {
    best_kappa <- NULL
    cv_mse <- NULL
  }

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

  v = .get_level_reg(xs, y, args$wts, args$nvars, args$nobs, args$q, intercept,
                     partial_method, args$l2_penalty, ...)
  meta_opt <- .get_meta_opt(y, args$kappa, args$nvars, args$nobs,
                            args$var_names_excl, standardize, intercept,
                            standard_sd, standard_mean, v)

  beta_mat <- meta_opt$beta
  opt_par_mat <- meta_opt$opt_par

  beta_mat_full <- array(NA, c(length(args$var_names), NCOL(beta_mat)),
                         dimnames = list(args$var_names, colnames(beta_mat)))
  beta_mat_full[rownames(beta_mat), colnames(beta_mat)] <- beta_mat

  if (intercept) fitted <- cbind(1, x) %*% beta_mat else fitted <- x %*% beta_mat
  resid <- y - fitted

  nlevels <- dim(v$fit_mat)[2]
  fit_mat_adj <- v$fit_mat
  theta_mat <- apply(opt_par_mat, 2, cumsum)
  for (i in 1:(nlevels-1)) fit_mat_adj[,i] <- fit_mat_adj[,i] - fit_mat_adj[,i+1]
  explained_variance <- sapply(nlevels:1, function(i) {
    if (i == nlevels) {
      fit <- fit_mat_adj[,nlevels] %*% matrix(theta_mat[nlevels,], nrow=1)
    } else {
      fit <- fit_mat_adj[,i:nlevels] %*% theta_mat[i:nlevels,]
    }
    return(1 - colSums(sweep(fit, 1, y)^2) / sum((y - mean(y))^2))
  })
  explained_variance <- t(explained_variance)

  out <- list(
    call = match.call(),
    coefficients = beta_mat_full,
    kappa = args$kappa,
    best_kappa = best_kappa,
    cv_mse = cv_mse,
    fitted.values = fitted,
    residuals = resid,
    x = x,
    y = y,
    df = round(as.numeric(v$dof %*% opt_par_mat), 4),
    hgraph = list(cluster_object = v$clust, shrinkage_vector = opt_par_mat,
                  included_levels = v$included_levels,
                  explained_variance = explained_variance,
                  full_level_output = v),
    intercept = intercept
  )

  class(out) <- "cv.hfr"

  return(out)

}
