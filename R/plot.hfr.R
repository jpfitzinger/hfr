#' @name plot.hfr
#' @title Plot the dendrogram of an HFR model
#' @description Plots the dendrogram of a fitted \code{hfr} model. The heights of the
#' dendrogram are given by a shrinkage vector, with a maximum (unregularized) height
#' of \eqn{(p-1)}{(p - 1)}. Stronger shrinkage leads to a shallower hierarchy.
#'
#' @details The dendrogram is generated using hierarchical clustering and modified
#' so that the height differential between any two splits is the shrinkage weight of
#' the lower split (ranging between 0 and 1). With no shrinkage, all shrinkage weights
#' are equal to 1 and the dendrogam has a height of \eqn{(p-1)}{(p - 1)}. With shrinkage
#' the dendrogram has a height of \eqn{(nu - )}{(nu - 1)}.
#'
#' The leaf nodes are colored to indicate the coefficient sign, with the size indicating
#' the absolute magnitude.
#'
#' @param x Fitted 'hfr' model.
#' @param ... additional methods passed to \code{plot}.
#' @return A plotted dendrogram.
#' @author Johann Pfitzinger
#'
#' @examples
#' x = matrix(rnorm(100 * 20), 100, 20)
#' y = rnorm(100)
#' fit = hfr(x, y, nu = 0.5)
#' plot(fit)
#'
#' @export
#'
#' @seealso \code{hfr}, \code{predict} and \code{coef} methods

plot.hfr <- function(
  x,
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

  heights <- rep(0, length(included_levels))
  heights[rev(included_levels)] <- theta

  var_names <- names(x$coefficients)
  if (x$intercept) var_names <- var_names[-1]

  .draw_dendro(clust, coefs, heights, x$hgraph$explained_variance, var_names, x$df)

}
