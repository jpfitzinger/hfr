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
#' The average standard errors along the branch of each coefficient can be used
#' to highlight coefficients that are not statistically significant. When
#' \code{confidence_level > 0}, branches with a lower confidence are plotted
#' as dotted lines.
#'
#' A color bar on the right indicates the relative contribution of each level to the
#' coefficient of determination, with darker hues representing a larger contribution.
#'
#' @param x Fitted 'hfr' model.
#' @param show_details print model details on the plot.
#' @param confidence_level coefficients with a lower approximate statistical confidence are highlighted in the plot, see details. Default is \code{confidence_level=0}.
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
#' @seealso \code{\link{hfr}}, \code{\link{se.avg}}, \code{\link{predict}} and \code{\link{coef}} methods
#'
#' @importFrom stats pt

plot.hfr <- function(
  x,
  show_details = TRUE,
  confidence_level = 0,
  max_leaf_size = 3,
  ...
) {

  if (!inherits(x, "hfr"))
    stop("object must be of class 'hfr'")

  if (confidence_level > 1 || confidence_level < 0)
    stop("'confidence_level' must be a scalar between 0 and 1")

  if (any(is.na(x$coefficients)))
    warning("removing variables with 'NA' coefficients")

  clust <- x$hgraph$cluster_object
  phi <- x$hgraph$shrinkage_vector
  included_levels <- x$hgraph$included_levels
  coefs <- x$coefficients
  coefs <- coefs[!is.na(coefs)]
  if (x$intercept) coefs <- coefs[-1]

  aggr <- diag(length(phi))
  aggr[lower.tri(aggr)] <- 1
  theta <- rev(as.numeric(aggr %*% phi))

  dof <- rev(sapply(x$hgraph$full_level_output$S, nrow))
  dof <- c(dof[1], dof[-1] - dof[-length(dof)])

  heights <- rep(0, length(included_levels))
  heights[rev(included_levels)] <- theta * dof

  var_names <- names(coefs)

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

  se <- se.avg(x)[-1]
  se <- se[!is.na(se)]
  pvals <- stats::pt(abs(coefs / se), NROW(x$y) - x$df - 1, lower.tail = F)
  dashed <- names(pvals)[pvals > 1 - confidence_level]

  .draw_dendro(clust, coefs, heights, expl_variance, var_names, x$df,
               show_details, max_leaf_size, dashed)

}
