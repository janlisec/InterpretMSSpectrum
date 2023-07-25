testthat::test_that(
  desc = "check_param_correct_peak works", 
  code = {
    inp <- "Name, C6H12O6, 0"
    testthat::expect_equal(as.character(InterpretMSSpectrum:::check_correct_peak(x = inp)), inp)
    # return NULL in case that input formula contains other elements than allowed by isotopes parameter
    testthat::expect_null(InterpretMSSpectrum:::check_correct_peak(x = inp, matrix("C", ncol=1)))
    testthat::expect_warning(InterpretMSSpectrum:::check_correct_peak(x = inp, matrix("C", ncol=1), silent=FALSE))
    testthat::expect_warning(InterpretMSSpectrum:::check_correct_peak(x = "wrong, string", silent=FALSE))
    testthat::expect_warning(InterpretMSSpectrum:::check_correct_peak(x = "wrong, string, really", silent=FALSE))
  }
)
