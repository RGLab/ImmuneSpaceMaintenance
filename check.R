suppressPackageStartupMessages(library(ImmuneSpaceMaintenance))
suppressPackageStartupMessages(library(Rlabkey))
suppressPackageStartupMessages(library(data.table))


login <- Sys.getenv("ISR_login")
password <- Sys.getenv("ISR_pwd")
machine <- Sys.getenv("ISR_machine")
# For dev purposes
machine <- ifelse(machine == "",
                  "test.immunespace.org",
                  machine)
string <- paste(
  "machine", machine,
  "login", login,
  "password", password
)
labkey.netrc.file <- tempfile()
write(string, labkey.netrc.file)

useragent <- paste0(
  "R/", R.version$major, ".", R.version$minor,
  " (", Sys.info()["sysname"], " ", Sys.info()["machine"], ")",
  " Rlabkey/", packageVersion("Rlabkey"))

dump <- labkey.setCurlOptions(
  ssl_verifyhost = 2,
  sslversion = 1,
  netrc_file = labkey.netrc.file,
  useragent = useragent
)

labkey.url.base <- paste0("https://", machine)
labkey.url.path <- "/Studies/"

check_type <- Sys.getenv("CHECK")
file_type <- Sys.getenv("FILE")
batch <- Sys.getenv("BATCH")

if (check_type == "checkStudyCompliance") {
  res <- checkStudyCompliance()

  # Remove studies with known GEM issues that are not fixable at the moment
  badGE <- c("SDY74",
             "SDY314",
             "SDY520",
             "SDY789",
             "SDY820",
             "SDY903",
             "SDY1092",
             "SDY1361")
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
      "\n"
    )
  }

} else if (check_type == "checkRawFiles" & file_type != "") {
  res <- checkRawFiles(file_type = file_type,
                             mc.cores = 1, # > 1 cores generates errors
                             batch = batch)

  if (sum(!res$file_exists) > 0) {
    res <- res[!res$file_exists, c("study_accession", "file_info_name")]
    res[ , list(files_missing = .N), by = study_accession]
    print(res)
    stop()
  }

} else if (check_type == "checkPublicVsStudySchema"){
  res <- checkPublicVsStudySchema()
  if (any(lengths(res)) > 0){
    res <- res[ lengths(res) > 0 ]
    print(res)
    stop()
  }
}
