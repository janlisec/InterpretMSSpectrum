testthat::test_that(
  desc = "findMAIN works", 
  code = {
    esi_spectrum <- InterpretMSSpectrum::esi_spectrum
    fmr <- InterpretMSSpectrum::findMAIN(esi_spectrum)
    
    # returns a object of class 'findMAIN'
    testthat::expect_true(inherits(fmr, "findMAIN"))
    
    # result should have 96 hypothesis tested (and be of length 96 therefore)
    testthat::expect_length(fmr, 96)
    
    # compute summary of fmr
    fmr_sum <- summary(fmr)
    
    # summary has 10 columns
    testthat::expect_equal(ncol(fmr_sum), 10)
    
    # summary has results for 3 main adduct hypotheses
    testthat::expect_true(all(fmr_sum[,"adducthyp"] %in% c("[M+H]+","[M+Na]+","[M+K]+")))
    
    # compute an alternative for only 4 peaks out of spec using only one adduct hypothesis
    fmr_sum2 <- summary(InterpretMSSpectrum::findMAIN(esi_spectrum[6:9,], adducthyp = "[M+H]+"))
    
    # summary2 has 2 rows
    testthat::expect_equal(nrow(fmr_sum2), 2)
    
    # summary2 has only M+H adduct hypotheses
    testthat::expect_true(all(fmr_sum2[,"adducthyp"] == c("[M+H]+")))
    
    # set up a spectrum containing a double charged peak
    spec <- data.frame(mz = c(372.1894, 372.6907, 373.1931, 380.1234), int = c(100, 40, 8, 2))
    
    # allow a double charged adduct hypothesis (not standard)
    fmr <- InterpretMSSpectrum::findMAIN(spec, adducthyp = c("[M+H]+", "[M+2H]2+"))
    fmr_sum3 <- summary(fmr)
    
    # check that correct adduct hyp is ranked first
    testthat::expect_true(fmr_sum3[1,2]=="[M+2H]2+")
    
    # add the correct M+H to this spectrum as a minor peak
    spec <- rbind(spec, c(742.3648+1.007, 10))
    fmr_sum4 <- summary(InterpretMSSpectrum::findMAIN(spec, adducthyp = c("[M+H]+", "[M+2H]2+")))
    
    # summary4 has 5 rows because one redundant hypothesis is removed
    testthat::expect_equal(nrow(fmr_sum4), 5)
    
    # compare specific hypotheses manually
    # get correct result and check total_score
    fmr_sum5 <- summary(InterpretMSSpectrum::findMAIN(spec, adductmz = 743.3718, adducthyp = "[M+H]+"))
    testthat::expect_equal(fmr_sum5[,10], 0.95)
    # enforce wrong result and check total_score
    fmr_sum6 <- summary(InterpretMSSpectrum::findMAIN(spec, adductmz = 743.3718, adducthyp = "[M+2H]2+"))
    testthat::expect_equal(fmr_sum6[,10], 0.51)
  }
)
