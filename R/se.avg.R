#' @name se.avg
#' @title Calculate approximate standard errors for a fitted HFR model
#' @description This function computes the weighted average standard errors across
#' levels using Burnham & Anderson (2004).
#'
#' @details The HFR computes linear regressions over several levels of an estimated
#' hierarchy. By averaging the standard errors across hierarchical levels, an
#' indication can be obtained about the average significance of the variables.
#'
#' Standard errors are understated, since the uncertainty in the hierarchy estimation
#' is not reflected.
#'
#' @param object Fitted \code{hfr} model.
#' @return A vector of standard errors.
#' @author Johann Pfitzinger
#' @references
#' Pfitzinger, J. (2022).
#' Cluster Regularization via a Hierarchical Feature Regression.
#' arXiv 2107.04831[statML]
#'
#' Burnham, K. P. and Anderson, D. R. (2004).
#' Multimodel inference - understanding AIC and BIC in model selection.
#' Sociological Methods & Research 33(2): 261-304.
#'
#' @examples
#' x = matrix(rnorm(100 * 20), 100, 20)
#' y = rnorm(100)
#' fit = hfr(x, y, kappa = 0.5)
#' se.avg(fit)
#'
#' @export
#'
#' @seealso \code{\link{hfr}} method
#'
#' @importFrom stats sd


se.avg <- function(
  object
) {

  if (!inherits(object, "hfr"))
    stop("object must be of class 'hfr'")

  standardize <- object$call$standardize
  intercept <- object$call$intercept

  if (is.null(standardize)) standardize <- TRUE else standardize <- as.logical(as.character(standardize))
  if (is.null(intercept)) intercept <- TRUE else intercept <- as.logical(as.character(intercept))

  nlevels <- length(object$hgraph$shrinkage_vector)

  coef_mat <- object$hgraph$full_level_output$coef_mat
  phi <- object$hgraph$shrinkage_vector

  beta_bar <- drop(coef_mat %*% phi)
  beta_sqdiff <- sweep(coef_mat, 1, beta_bar)^2

  lvl_reg <- object$hgraph$full_level_output$mod_list
  S <- object$hgraph$full_level_output$S

  if (intercept) {
    lvl_reg_var <- sapply(1:(nlevels-1), function(i) c(lvl_reg[[i]]$stderr[1], abs(t(S[[i]])) %*% lvl_reg[[i]]$stderr[-1]))^2
    lvl_reg_var <- cbind(lvl_reg_var, c(lvl_reg[[nlevels]]$stderr^2, rep(0, nrow(lvl_reg_var)-1)))
  } else {
    lvl_reg_var <- sapply(1:nlevels, function(i) abs(t(S[[i]])) %*% lvl_reg[[i]]$stderr)^2
  }


  stderr <- drop(sqrt(pmax(0, (lvl_reg_var + beta_sqdiff) %*% phi)))
  if (intercept) stderr <- stderr[-1]

  if (standardize) {
    standard_sd <- apply(object$x, 2, stats::sd)
    stderr <- stderr / standard_sd
  }

  coefs <- object$coefficients
  if (intercept) coefs <- coefs[-1]
  stderr_full <- stats::setNames(rep(NA, length(coefs)), names(coefs))
  stderr_full[!is.na(coefs)] <- stderr

  return(stderr_full)

}
