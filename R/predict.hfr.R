#' @name predict.hfr
#' @title Model predictions
#' @description Predict values using a fitted \code{hfr} model
#'
#' @details Predictions are made by multiplying the \code{newdata} object with the estimated coefficients.
#'
#' @param object Fitted 'hfr' model.
#' @param newdata Matrix or data.frame of new values for \code{x} at which predictions are to be made.
#' @param ... additional methods passed to \code{predict}.
#' @return A vector of predicted values.
#' @author Johann Pfitzinger
#'
#' @examples
#' x = matrix(rnorm(100 * 20), 100, 20)
#' y = rnorm(100)
#' fit = hfr(x, y, kappa = 0.5)
#' predict(fit)
#'
#' @export
#'
#' @seealso \code{\link{hfr}}, \code{\link{cv.hfr}} and \code{\link{coef}} methods
#'
#' @importFrom stats coef
#' @importFrom stats fitted

predict.hfr <- function(
  object,
  newdata = NULL,
  ...
  ) {

  if (!class(object) %in% c('hfr'))
    stop("object must be of class 'hfr'")

  if (is.null(newdata))
    return(stats::fitted(object))

  if (is.null(nobs <- nrow(newdata)))
    stop("'newdata' must be a matrix")
  if (nobs == 0L)
    stop("0 (non-NA) cases")
  nvars <- ncol(newdata)
  if (nvars != length(stats::coef(object)) - object$intercept)
    stop("incorrect number of columns in 'newdata'")

  if (any(is.na(newdata)))
    stop("'NA' values in 'newdata'")

  newdata <- data.matrix(newdata)

  if (object$intercept) {
    newdata <- cbind(1, newdata)
  }

  pred <- as.numeric(newdata %*% stats::coef(object))

  return(pred)

}
