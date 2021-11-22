#' @name predict.cv.hfr
#' @title Model Predictions
#' @description Predict values using a fitted Hierarchical Feature Regression cross-validation object
#'
#' @details The chosen hyperparameter value to use for predictions can be passed to
#' the \code{penalty} or \code{factors} argument.
#'
#' @param object Fitted 'cv.hfr' model.
#' @param newdata Matrix or data.frame of new values for \code{x} at which predictions are to be made.
#' @param penalty The optimal penalty used for prediction. Only when \code{object} is of type 'cv.hfr'.
#' @param factors The optimal factors used for prediction. Only when \code{object} is of type 'cv.hfr'.
#' @param ... additional methods passed to \code{predict}.
#' @return A vector of predicted values.
#' @author Johann Pfitzinger
#'
#' @examples
#' x = matrix(rnorm(100 * 20), 100, 20)
#' y = rnorm(100)
#' fit = cv.hfr(x, y, factors_grid = seq(0, 1, by = 0.1))
#' predict(fit, factors = 0.1)
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
  penalty = NULL,
  factors = NULL,
  ...
) {

  if (!class(object) %in% c('cv.hfr'))
    stop("object must be of class 'cv.hfr'")
  if (is.null(penalty) && is.null(factors))
    stop("must provide one of 'penalty' or 'factors'")
  if (!is.null(penalty)) {
    if (is.null(object$penalty_grid))
      stop("no 'penalty_grid' in 'object'")
    if (!any(penalty==object$penalty_grid))
      stop("'penalty' must be in 'penalty_grid' of the object")
    return_ix <- which(penalty==object$penalty_grid)
  }
  if (!is.null(factors)) {
    if (is.null(object$factors_grid))
      stop("no 'factors_grid' in 'object'")
    if (!any(factors==object$factors_grid))
      stop("'factors' must be in 'factors_grid' of the object")
    return_ix <- which(factors==object$factors_grid)
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

  if (object$intercept) {
    newdata <- cbind(1, newdata)
  }

  coefs <- stats::coef(object)[, return_ix]
  pred <- as.numeric(newdata %*% coefs)

  return(pred)

}
