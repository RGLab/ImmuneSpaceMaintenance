# ImmuneSpace Maintenance Tasks

<!-- badges: start -->
| [Production](https://www.immunespace.org/) | [Test](https://test.immunespace.org/) |
|-----|-----|
| [![IS-Maintenance](https://github.com/RGLab/ImmuneSpaceMaintenance/workflows/IS-Maintenance/badge.svg?branch=main)](https://github.com/RGLab/ImmuneSpaceMaintenance/actions?query=branch:main) | [![IS-Maintenance](https://github.com/RGLab/ImmuneSpaceMaintenance/workflows/IS-Maintenance/badge.svg?branch=dev)](https://github.com/RGLab/ImmuneSpaceMaintenance/actions?query=branch:dev) |
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
