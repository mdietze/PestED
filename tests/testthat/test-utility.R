context('utility functions')

test_that("format_inputs", {

  testthat::expect_error(format_inputs('fake.nc'))
  x <- format_inputs('met_data.nc')
  testthat::expect_true(is.data.frame(x))

})
