# hfr

`hfr` is an R package that implements the Hierarchical Feature Regression: a regularized regression estimator based on semi-supervised hierarchical graphs. The algorithm shrinks linear regression parameters towards group targets based on the similarity of the covariates.

## Usage

```
library(hfr)
x = matrix(rnorm(100 * 20), 100, 20)
y = rnorm(100)
fit = hfr(x, y, factors = 0.5)
coef(fit)
```

## Installation

```
devtools::install_github("https://github.com/jpfitzinger/hfr")
```

## References

Pfitzinger, J. (2021).
Cluster Regularization via a Hierarchical Feature Regression.
_arXiv 2107.04831[statML]_
