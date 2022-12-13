# ImmuneSpace Maintenance Tasks

<!-- badges: start -->
| [Production](https://datatools.immunespace.org/) | [Test](https://datatools-dev.immunespace.org/) |
|-----|-----|
| [![R-CMD-check](https://github.com/RGLab/ImmuneSpaceMaintenance/workflows/R-CMD-check/badge.svg?branch=main)](https://github.com/RGLab/ImmuneSpaceMaintenance/actions/workflows/R-CMD-check.yaml?query=branch:main) [![IS-Maintenance](https://github.com/RGLab/ImmuneSpaceMaintenance/workflows/IS-Maintenance/badge.svg?branch=main)](https://github.com/RGLab/ImmuneSpaceMaintenance/actions/workflows/IS-Maintenance.yaml?query=branch:main) | [![R-CMD-check](https://github.com/RGLab/ImmuneSpaceMaintenance/workflows/R-CMD-check/badge.svg?branch=dev)](https://github.com/RGLab/ImmuneSpaceMaintenance/actions/workflows/R-CMD-check.yaml?query=branch:dev) [![IS-Maintenance](https://github.com/RGLab/ImmuneSpaceMaintenance/workflows/IS-Maintenance/badge.svg?branch=dev)](https://github.com/RGLab/ImmuneSpaceMaintenance/actions/workflows/IS-Maintenance.yaml?query=branch:dev) |
<!-- badges: end -->

## Installation

``` r
# install.packages("remotes")
remotes::install_github("ImmuneSpaceMaintenance")
```

## Initialize an `ISM` connection

``` r
library(ImmuneSpaceMaintenance)
con <- ISM$new(study = "", onTest = TRUE)
```
