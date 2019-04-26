#' @include ISM.R
NULL



# PUBLIC -----------------------------------------------------------------------



# PRIVATE ----------------------------------------------------------------------

# Generate named list of files in either rawdata or analysis/exprs_matrices folders
ISM$set(
  which = "private",
  name = ".getGEFileNames",
  value = function(rawdata) {
    studies <- private$.getSdyVec()

    # check webdav folder for presence of rawdata
    file_list <- lapply(studies, FUN = function(sdy) {
      suffix <- ifelse(rawdata,
                       "/%40files/rawdata/gene_expression?method=JSON",
                       "/%40files/analysis/exprs_matrices?method=JSON"
      )

      dirLink <- paste0(
        self$config$labkey.url.base,
        "/_webdav/Studies/",
        sdy,
        suffix
      )
      files <- private$.listISFiles(dirLink)

      if (rawdata) {
        if (!is.null(files)) {
          files <- files[ grep("\\.(tsv|csv|cel|txt)$", files, ignore.case = T) ]
          files <- length(files) > 0
        }
      }

      return(files)
    })

    names(file_list) <- studies

    return(file_list)
  }
)


# Get vector of study folders
ISM$set(
  which = "private",
  name = ".getSdyVec",
  value = function() {
    studies <- labkey.getFolders(
      baseUrl = self$config$labkey.url.base,
      folderPath = "/Studies/"
    )[, 1]
    studies <- studies[grepl("SDY[0-9]+", studies)]

    studies
  }
)



# HELPER -----------------------------------------------------------------------
