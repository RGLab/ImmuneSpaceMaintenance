# NOTES:
# Bioconductor pkg 'ncdfFlow' runs into build error on linux os because
# compiler cannot find file.  Maybe -I flag is missing?
# In this case, need to use `osx` for operating system, but then
# `data.table` package fails being built from source due to different compiler
# error as osx doesn't use 'fopenmp'. Therefore need to pull binary for `data.table`
# from CRAN, which is only doable here as `r_binary_packages` is not usable with
# osx in travis.yml.

# Non-dependencies of package used here that must be installed
install.packages(c("devtools","testthat"))

# Installation of ISM and dependencies
devtools::install()
