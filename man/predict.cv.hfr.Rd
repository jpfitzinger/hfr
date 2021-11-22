% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict.cv.hfr.R
\name{predict.cv.hfr}
\alias{predict.cv.hfr}
\title{Model Predictions}
\usage{
\method{predict}{cv.hfr}(object, newdata = NULL, penalty = NULL, factors = NULL)
}
\arguments{
\item{object}{Fitted 'cv.hfr' model.}

\item{newdata}{Matrix or data.frame of new values for \code{x} at which predictions are to be made.}

\item{penalty}{The optimal penalty used for prediction. Only when \code{object} is of type 'cv.hfr'.}

\item{factors}{The optimal factors used for prediction. Only when \code{object} is of type 'cv.hfr'.}
}
\value{
A vector of predicted values.
}
\description{
Predict values using a fitted Hierarchical Feature Regression cross-validation object
}
\details{
The chosen hyperparameter value to use for predictions can be passed to
the \code{penalty} or \code{factors} argument.
}
\examples{
x = matrix(rnorm(100 * 20), 100, 20)
y = rnorm(100)
fit = cv.hfr(x, y, factors_grid = seq(0, 1, by = 0.1)
predict(fit, factors = 0.1)

}
\references{

}
\seealso{
\code{hfr}, \code{cv.hfr} and \code{coef} methods
}
\author{
Johann Pfitzinger
}