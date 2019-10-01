#' @include ISM.R
NULL



# PUBLIC -----------------------------------------------------------------------

# Remove gene expression matrices that do not correspond to a run currently on
# PROD/TEST in the query assay.ExpressionMatrix.matrix.Runs.
# NOTE: Important to change the labkey.url.base variable depending on PROD/TEST
# to ensure you are not deleting any incorrectly.
#' @importFrom httr DELETE HEAD http_type http_error
ISM$set(
  which = "public",
  name = "removeOrphanGEMs",
  value = function() {
    # Double check we are working at project level and on correct server!
    if (!private$.isProject()) {
      stop("Can only be run at project level")
    }
    chkBase <- readline(
      prompt = paste0(
        "You are working on ",
        self$config$labkey.url.base,
        ". Continue? [T / f] "
      )
    )
    if (!(chkBase %in% c("T", "t", ""))) {
      return("Operation Aborted.")
    }

    # get runs listed in the proper table
    runs <- data.table(
      labkey.selectRows(
        baseUrl = self$config$labkey.url.base,
        folderPath = self$config$labkey.url.path,
        schemaName = "assay.ExpressionMatrix.matrix",
        queryName = "Runs",
        showHidden = TRUE
      )
    )

    noRunPres <- private$.getNoRunPres(runs)

    # if files to-be-rm, confirm rm ok, attempt delete, check results and report
    if (length(noRunPres) == 0) {
      return("No orphans found.")
    } else {
      print(noRunPres)
    }
    ok2rm <- readline(prompt = "Ok to remove all files listed above? [Y / n] ")
    if (toupper(ok2rm) == "Y" | ok2rm == "") {
      for (i in seq_len(length(noRunPres))) {
        private$.curlDelete(
          noRunPres[i],
          names(noRunPres)[i]
        )
      }

      noRunPresPost <- private$.getNoRunPres(runs)
      if (length(noRunPresPost) == 0) {
        return("No orphans found after removal. Success!")
      } else {
        print("Problems Occurred. Remaining Files with No Runs Present")
        print(noRunPresPost)
        print("************")
      }
    } else {
      print("Operation aborted.")
    }
  }
)



# PRIVATE ----------------------------------------------------------------------

# Get no run pres (?)
ISM$set(
  which = "private",
  name = ".getNoRunPres",
  value = function(runs) {
    emFls <- private$.getGEFileNames(FALSE)
    emFls <- emFls[emFls != "NULL"]
    tmpNms <- rep(x = names(emFls), times = lengths(emFls))
    emFls <- unlist(emFls)
    names(emFls) <- tmpNms

    # get names from emFls for comparison
    emNms <- sapply(emFls, FUN = function(x) {
      return(strsplit(x, "\\.tsv")[[1]][1])
    })

    # check runs for file names - assume if not in runs$Name then ok to delete!
    noRunPres <- emNms[!(emNms %in% runs$Name)]
    noRunPres <- noRunPres[!duplicated(noRunPres)]

    noRunPres
  }
)


# Delete the files
ISM$set(
  which = "private",
  name = ".curlDelete",
  value = function(baseNm, sdy) {
    opts <- self$config$curlOptions
    opts$options$netrc <- 1L

    tsv <- paste0(
      self$config$labkey.url.base,
      "/_webdav/Studies/",
      sdy,
      "/%40files/analysis/exprs_matrices/",
      baseNm,
      ".tsv"
    )
    smry <- paste0(tsv, ".summary")

    tsvRes <- DELETE(url = tsv, config = opts)
    smryRes <- DELETE(url = smry, config = opts)

    list(tsv = tsvRes, summary = smryRes)
  }
)



# HELPER -----------------------------------------------------------------------
