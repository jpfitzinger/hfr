# hfr

`hfr` is an R package that implements a novel cluster-based regularized regression estimator: the **Hierarchical Feature Regression (HFR)**, which mobilizes insights from the domains of machine learning and graph theory to estimate robust parameters for a linear regression. The estimator constructs a supervised feature graph that decomposes parameters along its edges, adjusting first for common variation and successively incorporating idiosyncratic patterns into the fitting process. The graph structure has the effect of shrinking parameters towards group targets, where the extent of shrinkage is governed by a hyperparameter, and group compositions as well as shrinkage targets are determined endogenously. The hyperparameter `kappa` can take on values between `0` and `1` and represents the normalized effective model size. When `kappa = 1` the regression is unregularized resulting in OLS parameters. When `kappa < 1` the model is reduced to an effective model size smaller than the original covariate dimension.

## Usage

Fitting an `hfr` regression for the built-in `mtcars` data set:

```
library(hfr)
data("mtcars")
x <- mtcars[, -1]
y <- mtcars[, 1]
fit <- hfr(x, y, kappa = 0.75)
print(fit)
```

The HFR offers rich resources for the visual exploration of the underlying effect structure in the data. The package includes a custom dendrogram visualizing the optimal regression graph. See `?plot.hfr` for details:

```
plot(fit)
```

Optimal hyperparameter using k-fold cross-validation:

```
set.seed(123)
cv.fit <- cv.hfr(x, y, kappa_grid = seq(0, 1, by = 0.05))
cat("\nOptimal kappa: ", cv.fit$best_kappa)
plot(cv.fit)
```

The package provides functionality for the estimation of approximate standard errors:

```
fit <- hfr(x, y, kappa = cv.fit$best_kappa)
se.avg(fit)
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
