#' @name predict.cv.hfr
#' @title Model Predictions
#' @description Predict values using a fitted Hierarchical Feature Regression cross-validation object
#'
#' @details The chosen hyperparameter value to use for predictions can be passed to
#' the \code{penalty} or \code{factors} argument.
#'
#' @param object Fitted 'cv.hfr' model.
#' @param newdata Matrix or data.frame of new values for \code{x} at which predictions are to be made.
#' @param nu The optimal factors used for prediction.
#' @param ... additional methods passed to \code{predict}.
#' @return A vector of predicted values.
#' @author Johann Pfitzinger
#'
#' @examples
#' x = matrix(rnorm(100 * 20), 100, 20)
#' y = rnorm(100)
#' fit = cv.hfr(x, y, nu_grid = seq(0, 1, by = 0.1))
#' predict(fit, nu = 0.1)
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
  nu = NULL,
  ...
) {

  if (!class(object) %in% c('cv.hfr'))
    stop("object must be of class 'cv.hfr'")
  if (is.null(nu) && is.null(object$best_nu))
    stop("must provide 'nu'")
  if (is.null(nu) && !is.null(object$best_nu)) {
    nu <- object$best_nu
  }
  if (!is.null(nu)) {
    if (is.null(object$nu_grid))
      stop("no 'nu_grid' in 'object'")
    if (!any(nu==object$nu_grid))
      stop("'nu' must be in 'nu_grid' of the object")
    return_ix <- which(nu==object$nu_grid)
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
