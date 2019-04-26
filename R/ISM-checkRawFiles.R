#' @include ISM.R
NULL



# PUBLIC -----------------------------------------------------------------------

# Returns a list of data frames where TRUE in file_exists column marks files
# that are accessible. This function is used for administrative purposes to
# check that the raw files are properly loaded and accessible to the users.
#' @importFrom rjson fromJSON
#' @importFrom parallel mclapply detectCores
#' @importFrom Rlabkey labkey.getFolders
ISM$set(
  which = "public",
  name = "checkRawFiles",
  value = function(what = c(
    "gene_expression_files",
    "fcs_sample_files",
    "fcs_control_files",
    "protocols",
    "gene_expression_matrices"
  ),
  mc.cores = 1) {
    ## HELPERS
    ..messageResults <- function(dataset, file_exists) {
      message(
        paste0(
          sum(file_exists),
          "/",
          length(file_exists),
          " ",
          dataset,
          " with valid links."
        )
      )
    }

    ..checkLinks <- function(dataset, folder) {
      res <- data.frame(
        file_info_name = NULL,
        study_accession = NULL,
        file_link = NULL,
        file_exists = NULL,
        stringsAsFactors = FALSE
      )

      if (dataset %in% self$availableDatasets$Name) {
        temp <- self$getDataset(dataset, original_view = TRUE)

        if (dataset == "fcs_control_files") {
          temp <- temp[, file_info_name := control_file]
          temp <- temp[, c("pid", "sid") := data.table::tstrsplit(participant_id, "\\.")]
          temp <- temp[, study_accession := paste0("SDY", sid)]
        }

        temp <- temp[!is.na(file_info_name)]
        temp <- unique(temp[, list(study_accession, file_info_name)])

        file_link <- paste0(
          self$config$labkey.url.base,
          "/_webdav/Studies/",
          temp$study_accession,
          "/%40files/rawdata/",
          folder,
          "/",
          sapply(temp$file_info_name, URLencode)
        )

        studies <- unique(temp$study_accession)
        folder_link <- paste0(
          self$config$labkey.url.base,
          "/_webdav/Studies/",
          studies,
          "/%40files/rawdata/",
          folder,
          "?method=JSON"
        )

        file_list <- unlist(
          mclapply(
            folder_link,
            private$.listISFiles,
            mc.cores = mc.cores
          )
        )

        file_exists <- temp$file_info_name %in% file_list

        res <- data.frame(
          file_info_name = temp$file_info_name,
          study_accession = temp$study_accession,
          file_link = file_link,
          file_exists = file_exists,
          stringsAsFactors = FALSE
        )

        ..messageResults(dataset, res$file_exists)
      }

      res
    }


    ## MAIN

    startTimeTotal <- Sys.time()

    ret <- list()
    what <- tolower(what)

    if ("gene_expression_files" %in% what) {
      startTime <- Sys.time()

      ret$gene_expression_files <- ..checkLinks(
        "gene_expression_files",
        "gene_expression"
      )

      endTime <- Sys.time()
      print(endTime - startTime)
    }

    if ("fcs_sample_files" %in% what) {
      startTime <- Sys.time()

      ret$fcs_sample_files <- ..checkLinks(
        "fcs_sample_files",
        "flow_cytometry"
      )

      endTime <- Sys.time()
      print(endTime - startTime)
    }

    if ("fcs_control_files" %in% what) {
      startTime <- Sys.time()

      ret$fcs_control_files <- ..checkLinks(
        "fcs_control_files",
        "flow_cytometry"
      )

      endTime <- Sys.time()
      print(endTime - startTime)
    }

    if ("protocols" %in% what) {
      startTime <- Sys.time()

      if (private$.isProject()) {
        folders_list <- labkey.getFolders(
          baseUrl = self$config$labkey.url.base,
          folderPath = "/Studies/"
        )
        folders <- folders_list[, 1]
        folders <- folders[!folders %in% c("SDY_template", "Studies")]
      } else {
        folders <- basename(self$config$labkey.url.path)
      }

      file_link <- paste0(
        self$config$labkey.url.base,
        "/_webdav/Studies/",
        folders,
        "/%40files/protocols/",
        folders,
        "_protocol.zip"
      )

      file_exists <- unlist(
        mclapply(
          file_link,
          private$.checkUrl,
          mc.cores = mc.cores
        )
      )

      ..messageResults("protocols", file_exists)

      ret$protocols <- data.frame(
        file_info_name = paste0(folders, "_protocol.zip"),
        study_accession = folders,
        file_link = file_link,
        file_exists = file_exists,
        stringsAsFactors = FALSE
      )

      endTime <- Sys.time()
      print(endTime - startTime)
    }

    if ("gene_expression_matrices" %in% what) {
      startTime <- Sys.time()

      suppressWarnings(
        mx <- ImmuneSpaceR:::.getLKtbl(
          con = self,
          schema = "assay.ExpressionMatrix.matrix",
          query = "Runs",
          colNameOpt = "rname"
        )
      )

      if (nrow(mx) > 0) {
        mxLinks <- paste0(
          self$config$labkey.url.base,
          "/_webdav/Studies/",
          mx$folder_name,
          "/@files/analysis/exprs_matrices/",
          mx$name,
          ".tsv"
        )

        file_exists <- unlist(
          mclapply(
            mxLinks,
            private$.checkUrl,
            mc.cores = mc.cores
          )
        )

        ..messageResults("gene_expression_matrices", file_exists)

        ret$gene_expression_matrices <- data.frame(
          file_info_name = paste0(mx$name, ".tsv"),
          study_accession = mx$folder_name,
          file_link = mxLinks,
          file_exists = file_exists,
          stringsAsFactors = FALSE
        )
      } else {
        ret$gene_expression_matrices <- data.frame(
          file_info_name = NULL,
          study_accession = NULL,
          file_link = NULL,
          file_exists = NULL,
          stringsAsFactors = FALSE
        )
      }

      endTime <- Sys.time()
      print(endTime - startTime)
    }

    endTimeTotal <- Sys.time()
    print("===========")
    print("TOTAL TIME:")
    print(endTimeTotal - startTimeTotal)

    ret
  }
)



# PRIVATE ----------------------------------------------------------------------

# Check if the url exists (is accessible)
ISM$set(
  which = "private",
  name = ".checkUrl",
  value = function(url) {
    opts <- self$config$curlOptions
    opts$options$netrc <- 1L

    res <- HEAD(url, config = opts)

    if (http_error(res)) {
      ret <- FALSE
    } else {
      if (http_type(res) == "application/json") {
        res <- GET(url, config = opts)
        cont <- httr::content(res)
        ret <- is.null(cont$exception)
      } else {
        ret <- TRUE
      }
    }

    ret
  }
)



# HELPER -----------------------------------------------------------------------
