suppressPackageStartupMessages(library(ImmuneSpaceR))
suppressPackageStartupMessages(library(ImmuneSpaceMaintenance))

labkey.netrc.file <- ImmuneSpaceR:::.get_env_netrc()
labkey.url.base <- ImmuneSpaceR:::.get_env_url()

con <- ISM$new("")
msg <- testthat::capture_message(res <- con$checkStudyCompliance())
res

testthat::expect_equal(length(res), 0, label = msg$message)
