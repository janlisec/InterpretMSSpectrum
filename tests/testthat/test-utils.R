testthat::test_that(
  desc = "renderSMILES works correctly",
  code = {
    testthat::skip_on_ci()
    testthat::skip_on_cran()
    vdiffr::expect_doppelganger(
      title = "renderSMILES kekulise=FALSE",
      fig = function() {
        pdf(NULL)
        smiles <- "OS(=O)(=O)c1ccc(cc1)C(CC(=O)O)CC(=O)O"
        plot.new()
        plot.window(xlim=c(0,200), ylim=c(0,100))
        InterpretMSSpectrum:::renderSMILES(smiles, kekulise=FALSE)
      }
    )
    vdiffr::expect_doppelganger(
      title = "renderSMILES kekulise=TRUE",
      fig = function() {
        pdf(NULL)
        smiles <- "OS(=O)(=O)c1ccc(cc1)C(CC(=O)O)CC(=O)O"
        plot.new()
        plot.window(xlim=c(0,200), ylim=c(0,100))
        InterpretMSSpectrum:::renderSMILES(smiles, kekulise=TRUE)
      }
    )
  }
)

testthat::test_that(
  desc = "square_subplot_coord returns 4 numeric values",
  code = {
    grDevices::png(filename = tempfile(), width = 800, height = 600)
    on.exit(grDevices::dev.off())
    
    graphics::plot(1:10, 1:10)
    coords <- InterpretMSSpectrum:::square_subplot_coord(5, 5, w = 0.2)
    
    testthat::expect_type(coords, "double")
    testthat::expect_length(coords, 4)
  }
)

testthat::test_that(
  desc = "square_subplot_coord keeps center at x for inner points",
  code = {
    grDevices::png(filename = tempfile(), width = 800, height = 600)
    on.exit(grDevices::dev.off())
    
    graphics::plot(0:10, 0:10)
    coords <- InterpretMSSpectrum:::square_subplot_coord(5, 5, w = 0.2)
    
    center_x <- (coords[1] + coords[2]) / 2
    testthat::expect_equal(center_x, 5, tolerance = 1e-6)
  }
)
