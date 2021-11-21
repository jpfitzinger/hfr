
#' @importFrom cluster agnes
#' @importFrom cluster diana
#' @importFrom RcppArmadillo fastLmPure

.get_level_reg <- function(x, y, nvars, nobs, cluster_method, q, intercept) {

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
  clust_height_diff <- clust$height - c(clust$height[1], clust$height[-length(clust$height)])
  select_splits <- which(clust_height_diff>=quantile(clust_height_diff, 1-q))
  select_splits <- unique(c(1, select_splits, nvars))

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
  fit_mat <- c()
  dof <- c()

  for (i in 1:n_reg) {

    if (intercept) {

      x_i <- cbind(1, X_list[[i]])
      mod <- fastLmPure(x_i, y)
      mod_coef <- mod$coefficients
      mod_coef[is.na(mod_coef)] <- 0
      coef_list[[i]] <- c(mod_coef[1], t(S[[i]]) %*% mod_coef[-1])
      dof <- c(dof, length(mod_coef[-1]))

    } else {

      x_i <- X_list[[i]]
      mod <- fastLmPure(x_i, y)
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
    clust = clust
  ))

}
