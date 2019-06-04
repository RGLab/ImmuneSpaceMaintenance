suppressPackageStartupMessages(library(ImmuneSpaceR))
suppressPackageStartupMessages(library(ImmuneSpaceMaintenance))

labkey.netrc.file <- ImmuneSpaceR:::.get_env_netrc()
labkey.url.base <- ImmuneSpaceR:::.get_env_url()

con <- ISM$new("")
msg <- testthat::capture_messages(res <- con$checkStudyCompliance())

if (length(res) > 0) {
  dt <- data.table::as.data.table(do.call(rbind, res))
  dt[, study := names(res)]
  data.table::setcolorder(dt, c("study", "modules"))
  dt

  stop(
    "\n",
    length(res), " studies are not compliant: ",
    paste(names(res), collapse = ", "),
    "\n",
    msg
  )
}
