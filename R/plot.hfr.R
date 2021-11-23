#' @name plot.hfr
#' @title Plot the dendrogram of an HFR model
#' @description Plots the dendrogram of a fitted \code{hfr} model. The heights of the
#' dendrogram are given by a shrinkage vector, with a maximum (unregularized) height
#' of \eqn{(p-1)}{(p - 1)}. Stronger shrinkage leads to a shallower hierarchy.
#'
#' @details The dendrogram is generated using hierarchical clustering and modified
#' so that the height differential between any two splits is the shrinkage weight of
#' the lower split (ranging between 0 and 1). With no shrinkage, all shrinkage weights
#' are equal to 1 and the dendrogam has a height of \eqn{(p-1)}{(p - 1)}.
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
#'
#' @importFrom stats as.dendrogram
#' @importFrom graphics plot

plot.hfr <- function(
  x,
  ...
) {

  if (class(x)!="hfr")
    stop("object must be of class 'hfr'")

  clust <- x$cluster_model$cluster_object
  phi <- x$cluster_model$shrinkage_vector
  included_levels <- x$cluster_model$included_levels

  aggr <- diag(length(phi))
  aggr[lower.tri(aggr)] <- 1
  theta <- rev(as.numeric(aggr %*% phi))

  heights <- rep(0, length(included_levels))
  heights[rev(included_levels)] <- theta

  clust$height <- cumsum(rev(heights[-1]))

  graphics::plot(as.dendrogram(clust))

}
