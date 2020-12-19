# ImmuneSpace Maintenance Tasks

<!-- badges: start -->
  [![R-CMD-check-plus](https://github.com/RGLab/ImmuneSpaceMaintenance/workflows/R-CMD-check-plus/badge.svg)](https://github.com/RGLab/ImmuneSpaceMaintenance/actions)
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
