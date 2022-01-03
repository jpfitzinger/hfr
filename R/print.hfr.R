#' @name print.hfr
#' @title Print an HFR model
#' @description Print summary statistics for a fitted \code{hfr} model
#'
#' @details The call that produced the object \code{x} is printed, following by a
#' \code{data.frame} of summary statistics, including the effective degrees of freedom
#' of the model, the R.squared and the regularization parameter.
#'
#' @param x Fitted \code{hfr} model.
#' @param ... additional methods passed to \code{print}.
#' @return Summary statistics of HFR model
#' @author Johann Pfitzinger
#'
#' @examples
#' x = matrix(rnorm(100 * 20), 100, 20)
#' y = rnorm(100)
#' fit = hfr(x, y, kappa = 0.5)
#' print(fit)
#'
#' @export
#'
#' @seealso \code{hfr}, \code{predict}, \code{plot} and \code{coef} methods
#'
#' @importFrom stats as.dendrogram

print.hfr <- function(x, ...) {

  cat("\nCall: ", deparse(x$call), "\n\n")

  R2 <- 1 - sum(x$residuals^2) / sum((x$y - mean(x$y))^2)
  out = data.frame(Df = x$df, R.squared = round(R2, 2), check.names = FALSE)
  out$kappa <- x$kappa

  print(out)

}
