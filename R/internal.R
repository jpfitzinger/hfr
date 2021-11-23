
#' @importFrom stats hclust
#' @importFrom stats cor
#' @importFrom stats quantile
#' @importFrom stats cutree
#' @importFrom RcppArmadillo fastLmPure
#' @importFrom quadprog solve.QP

.get_level_reg <- function(x, y, nvars, nobs, q, intercept, ...) {

  corr <- stats::cor(x)
  cory <- stats::cor(x, y)

  rsq <- function(i, j) {

    rsq <- (cory[i,] - cory[j,] * corr[i,j]) / sqrt((1 - cory[j,]^2) * (1 - corr[i,j]^2))
    return(rsq)

  }
  partr2 <- matrix(NA, nvars, nvars)
  for (i in 1:nvars) for (j in 1:nvars) if (i != j) partr2[i,j] <- rsq(i, j)

  # Create distance matrix
  distmat <- stats::dist(partr2)

  clust <- hclust(distmat, ...)

  clust$height <- sort(clust$height)
  clust_height_diff <- clust$height - c(clust$height[1], clust$height[-length(clust$height)])
  select_splits <- which(clust_height_diff>=stats::quantile(clust_height_diff, 1-q))

  all_cuts_unfiltered <- c(nvars:1)
  all_cuts <- all_cuts_unfiltered[select_splits]
  all_cuts <- unique(c(min(nvars, nobs-2), all_cuts, 1))
  S <- c()
  X_list <- c()

  for (i in 1:length(all_cuts)) {

    cut <- stats::cutree(clust, k = all_cuts[i])
    max_cut <- max(cut)
    cut_fx <- function(row) as.numeric(row == cut)
    S_i <- lapply(1:max_cut, cut_fx)
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
  fit_mat <- c()
  dof <- c()

  for (i in 1:n_reg) {

    if (intercept) {

      x_i <- cbind(1, X_list[[i]])
      mod <- RcppArmadillo::fastLmPure(x_i, y)
      mod_coef <- mod$coefficients
      mod_coef[is.na(mod_coef)] <- 0
      coef_list[[i]] <- c(mod_coef[1], t(S[[i]]) %*% mod_coef[-1])
      dof <- c(dof, length(mod_coef[-1]))

    } else {

      x_i <- X_list[[i]]
      mod <- RcppArmadillo::fastLmPure(x_i, y)
      mod_coef <- mod$coefficients
      mod_coef[is.na(mod_coef)] <- 0
      coef_list[[i]] <- t(S[[i]]) %*% mod_coef
      dof <- c(dof, length(mod_coef))

    }

    fit_mat <- cbind(fit_mat, x_i %*% mod$coefficients)
    mod_list[[i]] <- mod

  }

  coef_mat <- sapply(coef_list, function(x) x)

  return(list(
    coef_mat = coef_mat,
    mod_list = mod_list,
    fit_mat = fit_mat,
    dof = dof,
    S = S,
    clust = clust,
    included_levels = all_cuts_unfiltered %in% all_cuts
  ))

}

.get_meta_opt <- function(y, nu, nvars, nobs, var_names, standardize, intercept, standard_sd, standard_mean, v) {

  grid_size <- length(nu)
  Dmat <- crossprod(v$fit_mat) / nobs
  diag(Dmat) <- diag(Dmat) + 1e-8
  dvec <- (t(v$fit_mat) %*% y) / nobs
  Amat <- diag(length(v$dof))
  # Amat[upper.tri(Amat)] <- 1
  Amat <- cbind(Amat, -Amat)
  bvec <- c(rep(0, length(v$dof)), -rep(1, length(v$dof)))

  beta_mat <- c()
  opt_par_mat <- c()
  for (i in 1:grid_size) {

    dof_constraint <- 1e-4 + nu[i] * (min(nvars, nobs-2)-1e-4-1e-8)
    opt <- quadprog::solve.QP(Dmat = Dmat,
                              dvec = dvec,
                              Amat = cbind(v$dof, Amat),
                              bvec = c(dof_constraint, bvec),
                              meq = 1)

    opt_par <- opt$solution

    beta <- rowSums(t(t(v$coef_mat) * opt_par))
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
    opt_par_mat <- cbind(opt_par_mat, opt_par)

  }

  colnames(beta_mat) <- nu

  return(list(beta = beta_mat, opt_par = opt_par_mat))

}
