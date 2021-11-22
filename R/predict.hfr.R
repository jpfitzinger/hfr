#' @name predict.hfr
#' @title Model Predictions
#' @description Predict values using a fitted Hierarchical Feature Regression
#'
#' @details Predictions are made by multiplying the newdata with the estimated coefficients.
#'
#' @param object Fitted 'hfr' or 'cv.hfr' model.
#' @param newdata Matrix or data.frame of new values for \code{x} at which predictions are to be made.
#' @return A vector of predicted values.
#' @author Johann Pfitzinger
#' @references
#'
#' @examples
#' x = matrix(rnorm(100 * 20), 100, 20)
#' y = rnorm(100)
#' fit = hfr(x, y, factors = 0.5)
#' predict(fit)
#'
#' @export
#'
#' @seealso \code{hfr}, \code{cv.hfr} and \code{coef} methods

predict.hfr <- function(
  object,
  newdata = NULL,
  penalty = NULL,
  factors = NULL
  ) {

  if (!class(object) %in% c('hfr'))
    stop("object must be of class 'hfr'")

  if (is.null(newdata))
    return(fitted(object))

  if (is.null(nobs <- nrow(newdata)))
    stop("'newdata' must be a matrix")
  if (nobs == 0L)
    stop("0 (non-NA) cases")
  nvars <- ncol(newdata)
  if (nvars != length(coef(object)) - object$intercept)
    stop("incorrect number of columns in 'newdata'")

  if (any(is.na(newdata)))
    stop("'NA' values in 'newdata'")

  if (object$intercept) {
    newdata <- cbind(1, newdata)
  }

  pred <- as.numeric(newdata %*% coef(object))

  return(pred)

}
