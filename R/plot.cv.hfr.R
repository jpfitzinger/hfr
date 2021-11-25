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
#' @importFrom graphics plot abline par
#' @importFrom dendextend set

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

  clust <- x$hgraph$cluster_object
  phi <- x$hgraph$shrinkage_vector[, return_ix]
  included_levels <- x$hgraph$included_levels
  coefs <- x$coefficients[, return_ix]
  if (x$intercept) coefs <- coefs[-1]
  coefs_sizes <- abs(coefs)/max(abs(coefs))
  coefs_sizes <- coefs_sizes * 3
  coefs_col <- ifelse(sign(coefs)>=0, "forestgreen", "firebrick")

  aggr <- diag(length(phi))
  aggr[lower.tri(aggr)] <- 1
  theta <- rev(as.numeric(aggr %*% phi))

  heights <- rep(0, length(included_levels))
  heights[rev(included_levels)] <- theta

  dend_heights <- cumsum(rev(heights[-1]))
  clust$height <- dend_heights
  var_names <- rownames(x$coefficients)
  if (x$intercept) var_names <- var_names[-1]

  dend <- stats::as.dendrogram(clust)
  dend <- dendextend::set(dend, "labels", var_names[clust$order])
  dend <- dendextend::set(dend, "leaves_pch", 15)
  dend <- dendextend::set(dend, "leaves_cex", coefs_sizes[clust$order])
  dend <- dendextend::set(dend, "leaves_col", coefs_col[clust$order])

  plot(x = rep(1, length(dend_heights)), y = dend_heights, type = "n", axes=F, xlab=NA, ylab=NA)
  for (i in dend_heights[dend_heights > 1e-4]) graphics::abline(h = i, col="lightgrey", lwd=1, lty = "dashed")
  graphics::par(new=TRUE)
  graphics::plot(as.dendrogram(dend))
  graphics::mtext(sprintf("Effective df: %.1f", x$df[return_ix]), side=3, line=1, at=0, col="black", las=1)

}
