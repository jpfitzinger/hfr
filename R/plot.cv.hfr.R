#' @name plot.cv.hfr
#' @title Plot the dendrogram of a cv.hfr object
#' @description Plots the dendrogram of a fitted \code{cv.hfr} model. The heights of the
#' dendrogram are given by a shrinkage vector, with a maximum (unregularized) height
#' of \eqn{(p-1)}{(p - 1)}. Stronger shrinkage leads to a shallower hierarchy.
#'
#' @details The dendrogram is generated using hierarchical clustering and modified
#' so that the height differential between any two splits is the shrinkage weight of
#' the lower split (ranging between 0 and 1). With no shrinkage, all shrinkage weights
#' are equal to 1 and the dendrogam has a height of \eqn{(p-1)}{(p - 1)}.
#'
#' @param x Fitted 'cv.hfr' model.
#' @param nu The optimal factors used for plotting.
#' @param ... additional methods passed to \code{plot}.
#' @return A plotted dendrogram.
#' @author Johann Pfitzinger
#'
#' @examples
#' x = matrix(rnorm(100 * 20), 100, 20)
#' y = rnorm(100)
#' fit = cv.hfr(x, y, nu_grid = seq(0, 1, by = 0.1))
#' plot(fit, nu = 0.5)
#'
#' @export
#'
#' @seealso \code{hfr}, \code{predict} and \code{coef} methods
#'
#' @importFrom stats as.dendrogram
#' @importFrom graphics plot

plot.cv.hfr <- function(
  x,
  nu = NULL,
  ...
) {

  if (!class(x) %in% c('cv.hfr'))
    stop("object must be of class 'cv.hfr'")
  if (is.null(nu) && is.null(x$best_nu))
    stop("must provide 'nu'")
  if (is.null(nu) && !is.null(x$best_nu)) {
    nu <- x$best_nu
  }
  if (!is.null(nu)) {
    if (is.null(x$nu_grid))
      stop("no 'nu_grid' in 'object'")
    if (!any(nu==x$nu_grid))
      stop("'nu' must be in 'nu_grid' of the object")
    return_ix <- which(nu==x$nu_grid)
  }

  clust <- x$cluster_model$cluster_object
  phi <- x$cluster_model$shrinkage_vector[, return_ix]
  included_levels <- x$cluster_model$included_levels

  aggr <- diag(length(phi))
  aggr[lower.tri(aggr)] <- 1
  theta <- rev(as.numeric(aggr %*% phi))

  heights <- rep(0, length(included_levels))
  heights[rev(included_levels)] <- theta

  clust$height <- cumsum(rev(heights[-1]))

  graphics::plot(as.dendrogram(clust))

}
