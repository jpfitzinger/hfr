## hfr 0.5.0

- Note that this starts from version `hfr 0.5.0`.

## hfr 0.6.0

- The plotting functionality has been expanded to highlight statistically insignificant branches
- Option to compute partial correlations using the pairwise (default) approach, or shrinkage approach
- Optional 'weights' argument to implement a weighted version of HFR.

## hfr 0.6.1

- Add an optional level-specific ridge_lambda parameter to penalize level-specific regressions
- Rename 'kappa_grid' argument to 'kappa' in 'cv.hfr' (for consistency across different methods)
- Minor bugfixes

## hfr 0.6.2

- Minor bugfixes in se.avg()

## hfr 0.7.0

- Introduce handling of linearly dependent columns and columns with zero standard deviation in 'x'
- Rename argument 'ridge_lambda' to 'l2_penalty'
- Updated documentation

## hfr 0.7.1

- fix handling of 'kappa' values < 1/k in the presence of intercept
