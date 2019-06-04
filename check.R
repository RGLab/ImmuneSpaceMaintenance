suppressPackageStartupMessages(library(ImmuneSpaceR))
suppressPackageStartupMessages(library(ImmuneSpaceMaintenance))

labkey.netrc.file <- ImmuneSpaceR:::.get_env_netrc()
labkey.url.base <- ImmuneSpaceR:::.get_env_url()

con <- ISM$new("")
res <- con$checkStudyCompliance()
res

testthat::expect_length(res, 0)
