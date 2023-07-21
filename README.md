
# InterpretMSSpectrum

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/InterpretMSSpectrum)](https://CRAN.R-project.org/package=InterpretMSSpectrum)
[![R-CMD-check](https://github.com/janlisec/InterpretMSSpectrum/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/janlisec/InterpretMSSpectrum/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The goal of **InterpretMSSpectrum** is to provides a set of R functions to annotate
mass spectra from Electrospray-Ionization and Atmospheric-Pressure-Chemical-Ionization 
derived data in positive and negative ionization mode.

## Installation

You can install the development version of InterpretMSSpectrum from 
[GitHub](https://github.com/) or
[CRAN](https://cran.r-project.org/package=InterpretMSSpectrum) with:

``` r
# devtools is required to install from GitHub
# install.packages("devtools")
devtools::install_github("janlisec/InterpretMSSpectrum")

# install the CRAN version als usual
install.packages("InterpretMSSpectrum")
```

## Example

In the simplest case **InterpretMSSpectrum** will provide an informed guess
for the potential sum formula of an arbitrary mass spectrum.

``` r
library(InterpretMSSpectrum)

# load APCI test data
apci_spectrum <- InterpretMSSpectrum::apci_spectrum

# find the most probable sum formula for the spectrum (will print to the console
# and open a new plot)
InterpretMSSpectrum(spec=apci_spectrum)
```

The function can be tweaked with numerous parameters to limit the results, speed
up calculations and more. The other high level function of the package allows
to predict the precursor of ESI spectra.

``` r
# load ESI test data
esi_spectrum <- InterpretMSSpectrum::esi_spectrum

# find the most likely precursor for the spectrum (will print to the console
# and open a new plot)
(fmr <- findMAIN(spec=esi_spectrum))
plot(fmr)
```

## Detailed documentation

You might read the publications on either
[APCI spectra processing](https://doi.org/10.1021/acs.analchem.6b02743) which
explains the idea of using the in source fragments for prediction of potential 
sum formulas or on
[ESI spectra processing](https://doi.org/10.1021/acs.analchem.6b02743) which 
explains the strategy to infer the correct precursor of ESI mass spectra.
