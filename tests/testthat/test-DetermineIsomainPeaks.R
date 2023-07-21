testthat::test_that(
  desc = "DetermineIsomainPeaks works", 
  code = {
    apci_spectrum <- InterpretMSSpectrum::apci_spectrum
    out1 <- InterpretMSSpectrum:::DetermineIsomainPeaks(spec=apci_spectrum, ionization="APCI")
    testthat::expect_length(out1, 5)
    testthat::expect_equal(out1, c(202.0893, 246.1339, 274.1288, 348.1477, 364.1789))
    
    esi_spectrum <- InterpretMSSpectrum::esi_spectrum
    out2 <- InterpretMSSpectrum:::DetermineIsomainPeaks(spec=esi_spectrum, ionization="ESI")
    testthat::expect_length(out2, 6)
    testthat::expect_equal(attr(out2, "pot_mh"), 263.05342)
    
  }
)
