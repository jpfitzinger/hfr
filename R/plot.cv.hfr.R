#' @name plot.cv.hfr
#' @title Plot the dendrogram of a cv.hfr object
#' @description Plots the dendrogram of a fitted \code{cv.hfr} model. The heights of the
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
#' @param x Fitted 'cv.hfr' model.
#' @param nu The optimal factors used for plotting.
#' @param show_details print model details on the plot.
#' @param max_leaf_size maximum size of the leaf nodes (default=3).
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

plot.cv.hfr <- function(
  x,
  nu = NULL,
  show_details = TRUE,
  max_leaf_size = 3,
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
    if (!any(round(nu, 6)==round(x$nu_grid, 6)))
      stop("'nu' must be in 'nu_grid' of the object")
    return_ix <- which(round(nu, 6)==round(x$nu_grid, 6))
  }

  clust <- x$hgraph$cluster_object
  phi <- x$hgraph$shrinkage_vector[, return_ix]
  included_levels <- x$hgraph$included_levels
  coefs <- x$coefficients[, return_ix]
  if (x$intercept) coefs <- coefs[-1]

  aggr <- diag(length(phi))
  aggr[lower.tri(aggr)] <- 1
  theta <- rev(as.numeric(aggr %*% phi))

  heights <- rep(0, length(included_levels))
  heights[rev(included_levels)] <- theta

  var_names <- rownames(x$coefficients)
  if (x$intercept) var_names <- var_names[-1]

  expl_variance <- rep(NA, length(included_levels))
  expl_variance[included_levels] <- x$hgraph$explained_variance[, return_ix]
  for (i in length(expl_variance):1) {
    if (is.na(expl_variance[i])) {
      if (i==length(expl_variance)) {
        expl_variance[i] <- 0
      } else {
        expl_variance[i] <- expl_variance[i+1]
      }
    }
  }

  .draw_dendro(clust, coefs, heights, expl_variance, var_names, x$df[return_ix], show_details, max_leaf_size)

}
