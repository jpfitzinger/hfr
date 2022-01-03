#' @name plot.hfr
#' @title Plot the dendrogram of an HFR model
#' @description Plots the dendrogram of a fitted \code{hfr} model. The heights of the
#' levels in the dendrogram are given by a shrinkage vector, with a maximum (unregularized)
#' overall graph height of \eqn{p}{p} (the number of covariates in the regression).
#' Stronger shrinkage leads to a shallower hierarchy.
#'
#' @details The dendrogram is generated using hierarchical clustering and modified
#' so that the height differential between any two splits is the shrinkage weight of
#' the lower split (ranging between \code{0} and \code{1}). With no shrinkage, all shrinkage weights
#' are equal to \code{1} and the dendrogram has a height of \eqn{p}{p}. With shrinkage
#' the dendrogram has a height of \eqn{(\kappa \times p)}{(\code{kappa} x p)}.
#'
#' The leaf nodes are colored to indicate the coefficient sign, with the size indicating
#' the absolute magnitude of the coefficients.
#'
#' A color bar on the right indicates the relative contribution of each level to the
#' coefficient of determination, with darker hues representing a larger contribution.
#'
#' @param x Fitted 'hfr' model.
#' @param show_details print model details on the plot.
#' @param max_leaf_size maximum size of the leaf nodes. Default is \code{max_leaf_size=3}.
#' @param ... additional methods passed to \code{plot}.
#' @return A plotted dendrogram.
#' @author Johann Pfitzinger
#'
#' @examples
#' x = matrix(rnorm(100 * 20), 100, 20)
#' y = rnorm(100)
#' fit = hfr(x, y, kappa = 0.5)
#' plot(fit)
#'
#' @export
#'
#' @seealso \code{hfr}, \code{predict} and \code{coef} methods

plot.hfr <- function(
  x,
  show_details = TRUE,
  max_leaf_size = 3,
  ...
) {

  if (class(x)!="hfr")
    stop("object must be of class 'hfr'")

  clust <- x$hgraph$cluster_object
  phi <- x$hgraph$shrinkage_vector
  included_levels <- x$hgraph$included_levels
  coefs <- x$coefficients
  if (x$intercept) coefs <- coefs[-1]

  aggr <- diag(length(phi))
  aggr[lower.tri(aggr)] <- 1
  theta <- rev(as.numeric(aggr %*% phi))

  dof <- rev(sapply(x$hgraph$full_level_output$S, nrow))
  dof <- c(dof[1], dof[-1] - dof[-length(dof)])

  heights <- rep(0, length(included_levels))
  heights[rev(included_levels)] <- theta * dof

  var_names <- names(x$coefficients)
  if (x$intercept) var_names <- var_names[-1]

  expl_variance <- rep(NA, length(included_levels))
  expl_variance[included_levels] <- x$hgraph$explained_variance
  for (i in length(expl_variance):1) {
    if (is.na(expl_variance[i])) {
      if (i==length(expl_variance)) {
        expl_variance[i] <- 0
      } else {
        expl_variance[i] <- expl_variance[i+1]
      }
    }
  }

  .draw_dendro(clust, coefs, heights, expl_variance, var_names, x$df, show_details, max_leaf_size)

}
