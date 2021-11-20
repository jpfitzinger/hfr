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
#' Hierarchical clustering is performed using the \code{cluster}-package. If
#' \code{cluster_method == 'DIANA'}, the function \code{cluster::diana} is used
#' to compute a cluster dendrogram, otherwise the function \code{cluster::agnes(., method = cluster_method)}
#' is used. Default is single-linkage agglomerative nesting.
#'
#' For high-dimensional problems, the hierarchy becomes very large. Setting \code{q} to a value below 1
#' reduces the number of levels used in the hierarchy. \code{q} represents a quantile-cutoff of amount of
#' information contributed by the levels. The default (\code{q = 1}) considers all levels.
#'
#' @param x Input matrix, of dimension \eqn{(N\times p)}{(N x p)}; each row is an observation vector.
#' @param y Response variable.
#' @param penalty A penalty on the effective degrees of freedom of the regression.
#' @param factors The target effective degrees of freedom of the regression.
#' @param q The quantile cut-off (in terms of information contributed) above which to consider levels in the hierarchy.
#' @param intercept Should intercept be fitted (default=TRUE).
#' @param standardize Logical flag for x variable standardization prior to fitting the model.
#' @param cluster_method hierarchical cluster algorithm used to construct an asset hierarchy.
#' @return An 'hfr' regression object.
#' @author Johann Pfitzinger
#' @references
#'
#' @examples
#' # Load returns of assets or portfolios
#' x = matrix(rnorm(100 * 20), 100, 20)
#' y = rnorm(100)
#' fit1 = glmnet(x, y)
#' print(fit1)
#'
#' @export
#'
#' @importFrom cluster agnes
#' @importFrom cluster diana
#' @importFrom quadprog solve.QP
#' @importFrom RcppArmadillo fastLmPure


hfr <- function(
  x,
  y,
  penalty = NULL,
  factors = NULL,
  q = NULL,
  intercept = TRUE,
  standardize = TRUE,
  cluster_method = c("DIANA", "single", "complete", "average", "ward")
  ) {

  cluster_method = match.arg(cluster_method)

  if (is.null(penalty) & is.null(factors)) {
    warning("Both 'penalty' and 'factors' are zero. Setting 'penalty = 0'")
    penalty <- 0
  }

  if (!is.null(factors)) {
    if (factors > 1 || factors < 0) {
      stop("'factors' must be between 0 and 1.")
    }
  }

  if (is.null(q)) {
    q <- min(nvars, sqrt(nobs)) / nvars
  }

  # Get feature names
  var_names <- colnames(x)
  if (is.null(var_names)) var_names <- paste("X", 1:ncol(x), sep = ".")
  if (intercept) var_names <- c("intercept", var_names)

  nobs = as.integer(dim(x)[1])
  nvars = as.integer(dim(x)[2])

  if (any(apply(x, 2, sd)==0))
    stop("Features can not have a standard deviation of zero.")

  if (standardize)
    x <- data.matrix(scale(x))

  corr <- cor(x)
  cory <- cor(x, y)

  rsq <- function(i, j) {

    rsq <- (cory[i,] - cory[j,] * corr[i,j]) / sqrt((1 - cory[j,]^2) * (1 - corr[i,j]^2))
    return(rsq)

  }
  partr2 <- matrix(NA, nvars, nvars)
  for (i in 1:nvars) for (j in 1:nvars) if (i != j) partr2[i,j] <- rsq(i, j)

  # Create distance matrix
  distmat <- dist(partr2)

  if (cluster_method == "DIANA") {
    clust <- diana(distmat)
  } else {
    clust <- agnes(distmat, method = cluster_method)
  }

  clust$height <- sort(clust$height)
  clust_height_diff <- clust$height[-1] - clust$height[-length(clust$height)]
  select_splits <- which(clust_height_diff>=quantile(clust_height_diff, 1-q))
  select_splits <- c(1, select_splits, nvars)

  all_cuts <- c(nvars:1)[select_splits]
  all_cuts <- all_cuts[all_cuts<=nobs-2]
  S <- c()
  X_list <- c()

  for (i in 1:length(all_cuts)) {

    cut <- cutree(clust, k = all_cuts[i])
    max_cut <- max(cut)
    cut_fx <- function(row) as.numeric(row == cut)
    S_i <- purrr::map(1:max_cut, cut_fx)
    S_i <- t(sapply(S_i, function(x) x))

    for (j in 1:nrow(S_i)) {

      if (sum(S_i[j,]==1) == 1) signs <- 1 else signs <- sign(colSums(corr[S_i[j,]==1, S_i[j,]==1]) - 1)
      signs[signs == 0] <- 1
      S_i[j,S_i[j,]==1] <- signs

    }

    S[[i]] <- S_i
    X_list[[i]] <- x %*% t(S[[i]])

  }

  n_reg <- length(S)
  coef_list <- c()
  mod_list <- c()

  for (i in 1:n_reg) {

    if (intercept) {

      mod <- fastLmPure(cbind(1, X_list[[i]]), y)
      mod_coef <- mod$coefficients
      mod_coef[is.na(mod_coef)] <- 0
      coef_list[[i]] <- c(mod_coef[1], t(S[[i]]) %*% mod_coef[-1])

    } else {

      mod <- fastLmPure(X_list[[i]], y)
      mod_coef <- mod$coefficients
      mod_coef[is.na(mod_coef)] <- 0
      coef_list[[i]] <- t(S[[i]]) %*% mod_coef

    }

    mod_list[[i]] <- mod

  }

  if (intercept) {
    fit_mat <- sapply(1:n_reg, function(i) cbind(1, x) %*% coef_list[[i]])
  } else {
    fit_mat <- sapply(1:n_reg, function(i) x %*% coef_list[[i]])
  }

  coef_mat <- sapply(coef_list, function(x) x)
  dof <- sapply(mod_list, function(x) length(coef(x)) - intercept)

  Dmat <- crossprod(fit_mat) * 2
  diag(Dmat) <- diag(Dmat) + 1e-8
  dvec <- 2*t(fit_mat) %*% y
  Amat <- diag(length(dof))
  Amat[upper.tri(Amat)] <- 1
  Amat <- cbind(Amat, -Amat)
  bvec <- c(rep(0, length(dof)), -rep(1, length(dof)))

  if (!is.null(penalty)) {
    penalty_term <- -dof * penalty * 100
    opt <- solve.QP(Dmat = Dmat,
                              dvec = dvec + penalty_term,
                              Amat = Amat,
                              bvec = bvec)
  } else {
    dof_constraint <- factors * nvars
    opt <- solve.QP(Dmat = Dmat,
                              dvec = dvec,
                              Amat = cbind(dof, Amat),
                              bvec = c(dof_constraint, bvec),
                              meq = 1)
  }

  opt.par <- opt$solution

  beta <- rowSums(t(t(coef_mat) * opt.par))
  names(beta) <- var_names

  if (intercept) fitted <- y - cbind(1, x) %*% beta else fitted <- y -  x %*% beta
  resid <- y - fitted

  out <- list(
    coefficients = beta,
    fitted = fitted,
    residuals = resid,
    theta = opt.par,
    dof = dof,
    fit_mat = fit_mat,
    coef_mat = coef_mat,
    S = S,
    clust = clust,
    BIC = nobs * log(mean(resid^2)) + 2 * log(nobs) * sum(dof * opt.par)
  )

  class(out) <- "hfr"

  return(out)

}
