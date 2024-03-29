#' @import ImmuneSpaceR data.table
#' @importFrom testthat capture_messages
#' @export
check <- function(check_type, file_type, batch = "") {
  con <- ISM$new("")

  if (check_type == "checkStudyCompliance") {
    msg <- testthat::capture_messages(
      res <- con$checkStudyCompliance()
    )

    # Remove studies with known GEM issues that are not fixable at the moment
    badGE <- c(
      "SDY74",
      "SDY314",
      "SDY789",
      "SDY820",
      "SDY903",
      "SDY1092",
      "SDY1361"
    )
    res <- res[!names(res) %in% badGE]

    if (length(res) > 0) {
      dt <- as.data.table(do.call(rbind, res))
      dt[, study := names(res)]
      setcolorder(dt, c("study", "modules"))
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
      res <- con$checkRawFiles(
        file_type = file_type,
        mc.cores = 1, # > 1 cores generates errors
        batch = batch
      )
    )

    if (sum(!res$file_exists) > 0) {
      res <- res[!res$file_exists, c("study_accession", "file_info_name")]
      res <- res[, list(files_missing = .N), by = study_accession]
      print(res)
      stop("Please download studies above!")
    }
    print(msg[1])
  } else if (check_type == "checkPublicVsStudySchema") {
    msg <- testthat::capture_messages(
      res <- con$checkPublicVsStudySchema()
    )
    if (any(lengths(res)) > 0) {
      res <- res[lengths(res) > 0]
      print(res)
      stop(msg)
    }
  }
}
