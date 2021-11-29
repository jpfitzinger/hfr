#' @name print.cv.hfr
#' @title Print a cv.hfr object
#' @description Print summary statistics for a fitted cv.hfr model
#'
#' @details The call that produced the object \code{x} is printed, following by a
#' data.frame of summary statistics, including the effective degrees of freedom
#' of the model, the R.squared and the regularization parameter.
#'
#' @param x Fitted 'cv.hfr' model.
#' @param ... additional methods passed to \code{print}.
#' @return Summary statistics of HFR model
#' @author Johann Pfitzinger
#'
#' @examples
#' x = matrix(rnorm(100 * 20), 100, 20)
#' y = rnorm(100)
#' fit = cv.hfr(x, y, nu_grid = seq(0, 1, by = 0.1))
#' print(fit)
#'
#' @export
#'
#' @seealso \code{hfr}, \code{predict}, \code{plot} and \code{coef} methods
#'
#' @importFrom stats as.dendrogram

print.cv.hfr <- function(x, ...) {

  cat("\nCall: ", deparse(x$call), "\n\n")

  cat("Best 'nu': ", x$best_nu, "\n\n")

  R2 <- c()
  for (i in 1:ncol(x$coefficients)) {
    R2 <- c(R2, 1 - sum(x$residuals[,i]^2) / sum(x$y^2))
  }
  out = data.frame(Df = x$df, R.squared = round(R2, 2), check.names = FALSE)
  out$nu <- x$nu_grid
  if (!is.null(x$cv_mse)) {
    out$MSE <- round(x$cv_mse, 2)
  }

  print(out)

}