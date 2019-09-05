suppressPackageStartupMessages(library(ImmuneSpaceR))
suppressPackageStartupMessages(library(ImmuneSpaceMaintenance))

labkey.netrc.file <- ImmuneSpaceR:::.get_env_netrc()
labkey.url.base <- ImmuneSpaceR:::.get_env_url()

check_type <- Sys.getenv("CHECK")
file_type <- Sys.getenv("FILE")

con <- ISM$new("")

if (check_type == "checkStudyCompliance") {
  msg <- testthat::capture_messages(
    res <- con$checkStudyCompliance()
  )

  # Remove studies with known GEM issues that are not fixable at the moment
  badGE <- c("SDY1092", "SDY74", "SDY820", "SDY1361", "SDY520", "SDY789")
  res <- res[ !names(res) %in% badGE ]

  if (length(res) > 0) {
    dt <- data.table::as.data.table(do.call(rbind, res))
    dt[, study := names(res)]
    data.table::setcolorder(dt, c("study", "modules"))
    print(dt[])

    stop(
      "\n",
      length(res), " studies are not compliant: ",
      paste(names(res), collapse = ", "),
      "\n",
      msg
    )
  }
} else if (check_type == "checkRawFiles" & file_type != "") {
  msg <- testthat::capture_messages(
    res <- con$checkRawFiles(file_type, mc.cores = parallel::detectCores())
  )
  res <- res[[file_type]]

  if (sum(!res$file_exists) > 0) {
    print(res[!res$file_exists, c("study_accession", "file_info_name")])
    stop(msg[1])
  }
}
