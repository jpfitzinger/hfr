# hfr

`hfr` is an R package that implements the Hierarchical Feature Regression: a regularized regression estimator based on semi-supervised hierarchical graphs.

## Installation

```
devtools::install_github("https://github.com/jpfitzinger/hfr")
```

## Usage

```
library(hfr)
x = matrix(rnorm(100 * 20), 100, 20)
y = rnorm(100)
fit = hfr(x, y, factors = 0.5)
coef(fit)
```

## References

Pfitzinger, J. (2021).
Cluster Regularization via a Hierarchical Feature Regression.
_arXiv_ [statML] 2107.04831
