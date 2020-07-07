
#' format_inputs
#'
#' @param nc_path Path to nc file
#' @return something
#' @export
format_inputs <- function(nc_path){

  assertthat::assert_that(file.exists(nc_path))

  # Import the data and format it
  nc <- ncdf4::nc_open(nc_path)

  PAR <- ncdf4::ncvar_get(nc, "PAR")
  for (i in which(PAR < -10)) {
    PAR[i] <- PAR[i - 1]
  } ## uber-naive gapfilling
  temp <- ncdf4::ncvar_get(nc, "TA")
  VPD <- ncdf4::ncvar_get(nc, "VPD")
  precip <- ncdf4::ncvar_get(nc, "PREC")
  time <- ncdf4::ncvar_get(nc, "DOY")
  ncdf4::nc_close(nc)

  data.frame(time = time, PAR = PAR, temp = temp, VPD = VPD, precip = precip)

}
