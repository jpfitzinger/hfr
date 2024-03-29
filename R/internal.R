
#' @importFrom stats hclust cor quantile cutree as.dendrogram as.dist dist
#' @importFrom quadprog solve.QP
#' @importFrom dendextend set get_nodes_xy branches_attr_by_labels
#' @importFrom RColorBrewer brewer.pal
#' @importFrom graphics abline segments points mtext rect plot
#' @importFrom corpcor pcor.shrink

.get_level_reg <- function(x, y, wts, nvars, nobs, q, intercept, partial_method, ridge_lambda, ...) {

  # Set default cluster method
  hclust_args <- match.call(expand.dots = F)$...
  if (!"method" %in% names(hclust_args)) {
    hclust_args$method <- "ward.D2"
  }

  corr <- stats::cor(x)
  cory <- stats::cor(x, y)

  if (partial_method == "pairwise") {

    partial_r_pairwise <- function(i, j) {

      coef <- (cory[i,] - cory[j,] * corr[i,j]) / sqrt((1 - cory[j,]^2) * (1 - corr[i,j]^2))
      return(coef)

    }
    partr <- matrix(NA, nvars, nvars)
    for (i in 1:nvars) for (j in 1:nvars) if (i != j) partr[i,j] <- partial_r_pairwise(i, j)

    distmat <- stats::dist(partr)

  } else if (partial_method == "shrinkage") {

    partial_r_shrinkage <- function(j) {
      out <- rep(NA, nvars)
      out[-j] <- drop(corpcor::pcor.shrink(cbind(y, x[,-j]), verbose = F)[1,-1])
      return(out)
    }

    partr <- sapply(1:ncol(x), partial_r_shrinkage)
    distmat <- stats::as.dist(sqrt((partr - t(partr))^2))

  }

  if (nvars > 2 & !any(is.na(distmat))) {
    clust <- do.call(hclust, append(list(d = distmat), hclust_args))
  } else {
    clust <- do.call(hclust, append(list(d = stats::dist(corr)), hclust_args))
  }

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
    cut_fx <- function(row) {
      ix <- row == cut
      if (sum(ix) == 1) {
        return(as.numeric(ix))
      } else {
        signs <- sign(drop(stats::cor(x, y)))
        signs[signs == 0] <- 1
        S_row <- as.numeric(ix) * signs
        return(S_row)
      }
    }
    S_i <- lapply(1:max_cut, cut_fx)
    S_i <- t(sapply(S_i, function(x) x))

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
      mod <- .custom_lm_ridge(x_i * wts, y * wts, lambda = ridge_lambda)
      mod_coef <- mod$coefficients
      mod_coef[is.na(mod_coef)] <- 0
      coef_list[[i]] <- c(mod_coef[1], t(S[[i]]) %*% mod_coef[-1])
      dof <- c(dof, length(mod_coef[-1]))

    } else {

      x_i <- X_list[[i]]
      mod <- .custom_lm_ridge(x_i * wts, y * wts, lambda = ridge_lambda)
      mod_coef <- mod$coefficients
      mod_coef[is.na(mod_coef)] <- 0
      coef_list[[i]] <- t(S[[i]]) %*% mod_coef
      dof <- c(dof, length(mod_coef))

    }

    fit_mat <- cbind(fit_mat, x_i %*% mod$coefficients)
    mod_list[[i]] <- mod

  }

  coef_mat <- sapply(coef_list, function(x) x)

  # add intercept level
  if (intercept) {
    coef_mat <- cbind(coef_mat, c(mean(y), rep(0, nvars)))
    mod_list[[nvars + 1]] <- list(coefficients = mean(y), stderr = sd(y)/sqrt(nobs))
    fit_mat <- cbind(fit_mat, mean(y))
    dof <- c(dof+1, 1)
  }

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

.get_meta_opt <- function(y, kappa, nvars, nobs, var_names, standardize, intercept, standard_sd, standard_mean, v) {

  grid_size <- length(kappa)
  Dmat <- crossprod(v$fit_mat) / nobs
  diag(Dmat) <- diag(Dmat) + 1e-8
  dvec <- (t(v$fit_mat) %*% y) / nobs
  Amat <- diag(length(v$dof))
  Amat <- cbind(Amat, -Amat)
  bvec <- c(rep(0, length(v$dof)), -rep(1, length(v$dof)))

  beta_mat <- c()
  opt_par_mat <- c()
  for (i in 1:grid_size) {

    dof_constraint <- intercept + kappa[i] * (min(nvars, nobs-2)-1e-8)
    for (eps in c(1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2)) {
      success <- tryCatch(
        {
          opt <- quadprog::solve.QP(Dmat = Dmat,
                             dvec = dvec,
                             Amat = cbind(v$dof, rep(-1, nvars+intercept), rep(1, nvars+intercept), Amat),
                             bvec = c(dof_constraint, -1-eps, 1-eps, bvec),
                             meq = 1)
          TRUE
        },
        error = function(cond) {
          FALSE
        }
      )
      if (success)
        break
    }
    if (!success)
      stop("optimization failed. check the input data.")


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

  colnames(beta_mat) <- paste0("s", 1:grid_size)

  return(list(beta = beta_mat, opt_par = opt_par_mat))

}

.draw_dendro <- function(clust, coefs, heights, explained_variance,
                         var_names, df, details, max_leaf_size,
                         dashed = NULL) {

  coefs_sizes <- abs(coefs)/max(abs(coefs))
  coefs_sizes <- coefs_sizes * max_leaf_size
  coefs_col <- ifelse(sign(coefs)>=0, "forestgreen", "firebrick")

  n <- length(heights)
  dend_heights <- cumsum(rev(heights))
  clust$height <- dend_heights[-n]

  dend <- stats::as.dendrogram(clust)
  dend <- dendextend::set(dend, "labels", var_names[clust$order])
  dend <- dendextend::set(dend, "leaves_pch", 15)
  dend <- dendextend::set(dend, "leaves_cex", coefs_sizes[clust$order])
  dend <- dendextend::set(dend, "leaves_col", coefs_col[clust$order])
  if (!is.null(dashed)) {
    dend <- dendextend::branches_attr_by_labels(dend, dashed, 3, "lty", "all")
  }

  pal <- RColorBrewer::brewer.pal(9, "Blues")
  pal <- c("#FFFFFF", pal)
  cols <- rev(explained_variance)
  cols <- pmax(c(cols[-n] - cols[-1], cols[n]), 0)
  cols <- round(sqrt((cols - min(cols)) / (max(cols) - min(cols))) * (length(pal)-1)+1)

  top_node <- dendextend::get_nodes_xy(dend)[1,]

  graphics::plot(stats::as.dendrogram(dend), ylim=c(0, max(dend_heights)),
                 panel.first=graphics::abline(h = dend_heights[dend_heights > 1e-4],
                                              col="lightgrey", lwd=1, lty = "dashed"))
  graphics::segments(x0 = top_node[1], y0 = top_node[2], y1 = dend_heights[n])
  graphics::points(x = top_node[1], y = dend_heights[n], pch = 15)
  graphics::rect(n*1.015, c(0, dend_heights[-n]), n*1.03, dend_heights, col = pal[cols], lwd=0.1)

  if (details) {
    graphics::mtext(sprintf("Effective df: %.1f; R-squared: %.2f", df, max(explained_variance)), side=3, line=1, at=0, adj=0, col="black", las=1)
  }

}

.custom_lm_ridge <- function(X, y, lambda) {

  p <- ncol(X)
  n <- nrow(X)
  XX <- crossprod(X)
  Xy <- crossprod(X, y)
  invXX <- solve(XX + lambda * n * diag(p))
  beta <- drop(invXX %*% Xy)
  if (n > p) {
    var_hat <- sum((y - X %*% beta)^2) / (n - p)
    std_err <- sqrt(diag(invXX * var_hat))
  } else {
    std_err <- rep(NA, length(beta))
  }
  return(list(coefficients = beta, stderr = std_err))

}
