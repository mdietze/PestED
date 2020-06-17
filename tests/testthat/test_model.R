
context('SEM default run')

test_that("iterate.SEM", {

  # Read in inputs
  x <- format_inputs('met_data.nc')
  suppressMessages({
    xx  <- iterate.SEM(pest = c(0, 0, 0, 0, 0), inputs = x, params = default_parameters, years = 1)
  })
  testthat::expect_true(is.data.frame(xx))
  testthat::expect_true(nrow(xx) > 1)


})
