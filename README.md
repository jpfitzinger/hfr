# hfr

`hfr` is an R package that implements the Hierarchical Feature Regression: a regularized regression estimator based on a supervised hierarchical graph. The algorithm identifies predictors with a similar effect on the response variable, and shrinks linear regression parameters towards group targets for these predictors. The strength of shrinkage is governed by a hyperparameter `kappa` that can take on values between 0 and 1 and represents the normalized effective model size. When `kappa = 1` the regression is unregularized resulting in OLS parameters. When `kappa < 1` the model is reduced to an effective model size of $1 + (p-1)\kappa$, where $p$ is the number of covariates.

## Usage

Fitting an `hfr` regression for the built-in `mtcars` data set:

```
library(hfr)
data("mtcars")
x <- mtcars[, -1]
y <- mtcars[, 1]
fit <- hfr(x, y, kappa = 0.75)
print(fit)
plot(fit)
```

Optimal hyperparameter using $k$-fold cross-validation:

```
set.seed(123)
cv.fit <- cv.hfr(x, y, kappa_grid = seq(0, 1, by = 0.05))
cat("\nOptimal kappa: ", cv.fit$best_kappa)
plot(cv.fit)
```

The package provides functionality of the estimation of approximate standard errors:

```
fit <- hfr(x, y, kappa = cv.fit$best_kappa)
se_avg(fit)
```

Finally, standard functions such as `coef`, `predict` and `print` can be used to interact with fitted `hfr` or `cv.hfr` objects.

## Installation

```
devtools::install_github("https://github.com/jpfitzinger/hfr")
```

## References

Pfitzinger, J. (2021).
Cluster Regularization via a Hierarchical Feature Regression.
_arXiv 2107.04831[statML]_
