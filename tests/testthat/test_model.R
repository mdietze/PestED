
context('SEM default run')

test_that("iterate.SEM", {

  # Read in inputs
  x <- format_inputs('met_data.nc')
  timestep <- 1800
  suppressMessages({
    xx  <- iterate.SEM(pest = c(0, 0, 0, 0, 0),
                       inputs = x,
                       params = default_parameters,
                       timestep = 1800,
                       years = 1)
  })
  expect_true(is.data.frame(xx))
  expect_true(nrow(xx) > 1)
})
