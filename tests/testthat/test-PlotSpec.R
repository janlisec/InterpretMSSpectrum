testthat::test_that(
  desc = "PlotSpec plot returns expected result",
  code = {
    inp <- InterpretMSSpectrum::apci_spectrum
    # avoid creating a Rplots.pdf in testthat folder
    pdf(NULL)
    vdiffr::expect_doppelganger(
      title = "PlotSpec_01",
      fig = function() InterpretMSSpectrum::PlotSpec(x = inp, ionization="APCI")
    )
  }
)
