.check_args <- function(x, y, intercept, weights, kappa, q, l2_penalty, is_cv) {
  nvars <- ncol(x)
  if (is.null(nobs <- nrow(x)))
    stop("'x' must be a matrix")
  if (nobs == 0L)
    stop("0 (non-NA) cases")
  ny <- NCOL(y)
  if (is.matrix(y) && ny == 1)
    y <- drop(y)
  if (ny > 1)
    stop("'y' must be a single response variable")
  if (NROW(y) != nobs)
    stop("incompatible dimensions")

  if (any(is.na(y)) || any(is.na(x)))
    stop("'NA' values in 'x' or 'y'")

  if (is_cv) {
    if (any(kappa > 1) || any(kappa < 0)) {
      stop("each 'kappa' must be between 0 and 1.")
    }
  } else {
    if (kappa > 1 || kappa < 0) {
      stop("'kappa' must be between 0 and 1")
    }
  }

  if (!is.null(weights)) {
    if (length(weights) != nobs)
      stop("'weights' must have same length as 'y'")
    if (any(is.na(weights)))
      stop("'NA' values in 'weights'")
    if (any(weights < 0))
      stop("'weights' can only contain positive numerical values")
    wts <- sqrt(weights)
  } else {
    wts <- rep(1, nobs)
  }

  # Get feature names
  var_names <- colnames(x)
  if (is.null(var_names)) var_names <- paste("V", 1:ncol(x), sep = "")
  if (intercept) var_names <- c("(Intercept)", var_names)

  # Convert 'x' to matrix
  x <- data.matrix(x)

  if (l2_penalty < 0)
    stop("'l2_penalty' must be a positive value")

  # Check for linearly dependent columns
  x_rank <- qr(crossprod(x))$rank
  if ((x_rank < nvars) & l2_penalty == 0) {
    warning("linearly dependent columns in 'x'. setting 'l2_penalty' to a small positive value")
    l2_penalty <- 1e-4
  }

  # Remove features with no standard deviation
  remove_columns <- which(apply(x, 2, stats::sd) == 0)
  if (length(remove_columns) > 0) {
    x <- x[, -remove_columns]
    var_names_excl <- var_names[-remove_columns-intercept]
    nvars <- ncol(x)
  } else {
    var_names_excl <- var_names
  }

  if (is.null(q)) {
    q <- 1
  }

  return(list(x = x, y = y, wts = wts, kappa = kappa, q = q, nobs = nobs, nvars = nvars,
              l2_penalty = l2_penalty, var_names = var_names, var_names_excl = var_names_excl))
}
