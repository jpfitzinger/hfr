#' @name predict.cv.hfr
#' @title Model predictions
#' @description Predict values using a fitted \code{cv.hfr} model
#'
#' @details Predictions are made by multiplying the \code{newdata} object with the estimated coefficients.
#' The chosen hyperparameter value to use for predictions can be passed to
#' the \code{kappa} argument.
#'
#' @param object Fitted 'cv.hfr' model.
#' @param newdata Matrix or data.frame of new values for \code{x} at which predictions are to be made.
#' @param kappa The hyperparameter used for prediction. If empty, the optimal value is used.
#' @param ... additional methods passed to \code{predict}.
#' @return A vector of predicted values.
#' @author Johann Pfitzinger
#'
#' @examples
#' x = matrix(rnorm(100 * 20), 100, 20)
#' y = rnorm(100)
#' fit = cv.hfr(x, y, kappa_grid = seq(0, 1, by = 0.1))
#' predict(fit, kappa = 0.1)
#'
#' @export
#'
#' @seealso \code{hfr}, \code{cv.hfr} and \code{coef} methods
#'
#' @importFrom stats fitted
#' @importFrom stats coef

predict.cv.hfr <- function(
  object,
  newdata = NULL,
  kappa = NULL,
  ...
) {

  if (!class(object) %in% c('cv.hfr'))
    stop("object must be of class 'cv.hfr'")
  if (is.null(kappa) && is.null(object$best_kappa))
    stop("must provide 'kappa'")
  if (is.null(kappa) && !is.null(object$best_kappa)) {
    kappa <- object$best_kappa
  }
  if (!is.null(kappa)) {
    if (is.null(object$kappa_grid))
      stop("no 'kappa_grid' in 'object'")
    if (!any(kappa==object$kappa_grid))
      stop("'kappa' must be in 'kappa_grid' of the object")
    return_ix <- which(kappa==object$kappa_grid)
  }

  if (is.null(newdata)) {
    return(stats::fitted(object)[,return_ix])
  }

  if (is.null(nobs <- nrow(newdata)))
    stop("'newdata' must be a matrix")
  if (nobs == 0L)
    stop("0 (non-NA) cases")
  nvars <- ncol(newdata)
  if (nvars != nrow(stats::coef(object)) - object$intercept)
    stop("incorrect number of columns in 'newdata'")

  if (any(is.na(newdata)))
    stop("'NA' values in 'newdata'")

  newdata <- data.matrix(newdata)

  if (object$intercept) {
    newdata <- cbind(1, newdata)
  }

  coefs <- stats::coef(object)[, return_ix]
  pred <- as.numeric(newdata %*% coefs)

  return(pred)

}
