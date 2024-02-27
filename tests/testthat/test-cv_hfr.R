test_that("cv_hfr", {
  data("Boston", package = "MASS")
  x <- model.matrix(medv ~ 0 + . + crim*tax, Boston)
  y <- Boston$medv
  expect_no_error(mod_1 <- cv.hfr(x, y))

  expect_no_error(coef(mod_1))
  expect_no_error(predict(mod_1, x))
  expect_no_error(fitted(mod_1))
  expect_error(se.avg(mod_1))
  expect_no_error(plot(mod_1))

  expect_false(any(is.na(coef(mod_1))))
  expect_false(any(is.na(predict(mod_1, x))))
  expect_false(any(is.na(fitted(mod_1))))
})

test_that("hfr-linear-dep", {
  data("Boston", package = "MASS")
  x <- model.matrix(medv ~ 0 + . + crim*tax, Boston)
  y <- Boston$medv
  x[,1] <- x[,2]
  expect_warning(mod_1 <- cv.hfr(x, y))

  expect_no_error(coef(mod_1))
  expect_no_error(predict(mod_1, x))
  expect_no_error(fitted(mod_1))
  expect_no_error(plot(mod_1))

  expect_false(any(is.na(coef(mod_1))))
  expect_false(any(is.na(predict(mod_1, x))))
  expect_false(any(is.na(fitted(mod_1))))
})

test_that("hfr-const-col", {
  data("Boston", package = "MASS")
  x <- model.matrix(medv ~ 0 + . + crim*tax, Boston)
  y <- Boston$medv
  x[,1] <- 1
  expect_no_error(mod_1 <- cv.hfr(x, y))

  expect_no_error(coef(mod_1))
  expect_no_error(predict(mod_1, x))
  expect_no_error(fitted(mod_1))
  expect_warning(plot(mod_1))

  expect_true(any(is.na(coef(mod_1))))
  expect_false(any(is.na(predict(mod_1, x))))
  expect_false(any(is.na(fitted(mod_1))))
})
