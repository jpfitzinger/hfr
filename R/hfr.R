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
#' When \code{kappa} is \code{1} the result is equivalent to the result from an ordinary
#' least squares regression (no shrinkage). Conversely, \code{kappa} set to \code{0}
#' represents maximum shrinkage.
#'
#' When \eqn{p > N}{p > N} \code{kappa} is a percentage of \eqn{(N - 2)}{(N - 2)}.
#'
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
#' When data exhibits multicollinearity it can be useful to include a penalty on the l2 norm in the level-specific regressions.
#' This can be achieved by setting the \code{l2_penalty} parameter.
#'
#' @param x Input matrix or data.frame, of dimension \eqn{(N\times p)}{(N x p)}; each row is an observation vector.
#' @param y Response variable.
#' @param weights an optional vector of weights to be used in the fitting process. Should be NULL or a numeric vector. If non-NULL, weighted least squares is used for the level-specific regressions.
#' @param kappa The target effective degrees of freedom of the regression as a percentage of \eqn{p}{p}.
#' @param q Thinning parameter representing the quantile cut-off (in terms of contributed variance) above which to consider levels in the hierarchy. This can used to reduce the number of levels in high-dimensional problems. Default is no thinning.
#' @param intercept Should intercept be fitted. Default is \code{intercept=TRUE}.
#' @param standardize Logical flag for x variable standardization prior to fitting the model. The coefficients are always returned on the original scale. Default is \code{standardize=TRUE}.
#' @param partial_method Indicate whether to use pairwise partial correlations, or shrinkage partial correlations.
#' @param l2_penalty Optional penalty for level-specific regressions (useful in high-dimensional case)
#' @param ...  Additional arguments passed to \code{hclust}.
#' @return An 'hfr' regression object.
#' @author Johann Pfitzinger
#' @references
#' Pfitzinger, Johann (2024). Cluster Regularization via a Hierarchical Feature Regression. _Journal of Econometrics and Statistics_ (in press). URL https://doi.org/10.1016/j.ecosta.2024.01.003.
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
#' @importFrom stats sd setNames


hfr <- function(
  x,
  y,
  weights = NULL,
  kappa = 1,
  q = NULL,
  intercept = TRUE,
  standardize = TRUE,
  partial_method = c("pairwise", "shrinkage"),
  l2_penalty = 0,
  ...
  ) {

  args = .check_args(x = x, y = y, intercept = intercept, kappa = kappa,
                     weights = weights, q = q, l2_penalty = l2_penalty, is_cv = FALSE)
  if (args$nvars == 0L) {
    return(list(coefficients = numeric(), residuals = y,
                fitted.values = 0 * y, dof = 0, clust = NULL,
                intercept = intercept))
  }
  partial_method = match.arg(partial_method)

  x <- args$x
  y <- args$y

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
  meta_opt <- .get_meta_opt(y, args$kappa, args$nvars, args$nobs, args$var_names_excl,
                            standardize, intercept, standard_sd, standard_mean, v)

  beta <- drop(meta_opt$beta)
  opt_par <- drop(meta_opt$opt_par)

  if (intercept) fitted <- as.numeric(cbind(1, x) %*% beta) else fitted <- as.numeric(x %*% beta)
  resid <- as.numeric(y - fitted)

  beta_full <- stats::setNames(rep(NA, length(args$var_names)), args$var_names)
  beta_full[args$var_names_excl] <- beta

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
    coefficients = beta_full,
    kappa = args$kappa,
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
