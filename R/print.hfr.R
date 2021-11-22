#' @name print.hfr
#' @title Print an HFR object
#' @description Print summary statistics for a fitted HFR model
#'
#' @details The call that produced the object \code{x} is printed, following by a
#' data.frame of summary statistics, including the effective degrees of freedom
#' of the model, the R.squared and the regularization parameter.
#'
#' @param object Fitted 'hfr' model.
#' @return Summary statistics of HFR model
#' @author Johann Pfitzinger
#' @references
#'
#' @examples
#' x = matrix(rnorm(100 * 20), 100, 20)
#' y = rnorm(100)
#' fit = hfr(x, y, factors = 0.5)
#' print(fit)
#'
#' @export
#'
#' @seealso \code{hfr}, \code{predict}, \code{plot} and \code{coef} methods
#'
#' @importFrom stats as.dendrogram

print.hfr <- function(object) {

  cat("\nCall: ", deparse(object$call), "\n\n")
  R2 <- 1 - sum(object$residuals^2) / sum(object$y^2)

  out = data.frame(Df = object$df, R.squared = round(R2, 2), check.names = FALSE)
  if (!is.null(object$penalty)) out$penalty <- object$penalty
  if (!is.null(object$factors)) out$factors <- object$factors
  print(out)

}
