#' @name plot.hfr
#' @title Plot the dendrogram of an HFR model
#' @description Plots the dendrogram of a fitted \code{hfr} model. The heights of the
#' dendrogram are given by a shrinkage vector, with a maximum (unregularized) height
#' of \eqn{(p-1)}{(p - 1)}. Stronger shrinkage leads to a shallower hierarchy.
#'
#' @details
#'
#' @param object Fitted 'hfr' model.
#' @return A plotted dendrogram.
#' @author Johann Pfitzinger
#' @references
#'
#' @examples
#' x = matrix(rnorm(100 * 20), 100, 20)
#' y = rnorm(100)
#' fit = hfr(x, y, factors = 0.5)
#' plot(fit)
#'
#' @export
#'
#' @seealso \code{hfr}, \code{predict} and \code{coef} methods
#'
#' @importFrom stats as.dendrogram

plot.hfr <- function(
  object
) {

  if (class(object)!="hfr")
    stop("object must be of class 'hfr'")

  clust <- object$cluster_model$cluster_object
  phi <- object$cluster_model$shrinkage_vector
  included_levels <- object$cluster_model$included_levels

  aggr <- diag(length(phi))
  aggr[lower.tri(aggr)] <- 1
  theta <- rev(as.numeric(aggr %*% phi))

  heights <- rep(0, length(included_levels))
  heights[rev(included_levels)] <- theta

  clust$height <- cumsum(rev(heights[-1]))

  plot(as.dendrogram(clust))

}
