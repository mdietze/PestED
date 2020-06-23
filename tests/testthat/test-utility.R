context('utility functions')

test_that("format_inputs", {

  expect_error(format_inputs('fake.nc'))
  x <- format_inputs('met_data.nc')
  expect_true(is.data.frame(x))

})
