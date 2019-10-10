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
  value = function(file_type, mc.cores, batch) {

    ## ------- HELPERS --------
    ..messageResults <- function(file_type, file_exists) {
      message(
        paste0(
          sum(file_exists),
          "/",
          length(file_exists),
          " ",
          file_type,
          " with valid links."
        )
      )
    }

    ..checkLinksRawFolder <- function(file_type, folder, batch) {
        temp <- self$getDataset(file_type, original_view = TRUE)

        if (file_type == "fcs_control_files") {
          temp <- temp[, file_info_name := control_file]
          temp <- temp[, c("pid", "sid") := data.table::tstrsplit(participant_id, "\\.")]
          temp <- temp[, study_accession := paste0("SDY", sid)]
        }

        temp <- temp[!is.na(file_info_name)]
        temp <- unique(temp[, list(study_accession, file_info_name)])

        # Batch system created as TravisCI has 50 min limit per job
        # but unlimited jobs. With increasing number of studies,
        # FCS file checking > 50 min at project level.
        # Assuming only 2 batches in this code
        if(batch != ""){
          mid <- dim(temp)[1]/2
          initSdy <- currSdy <- temp$study_accession[mid]
          while(initSdy == currSdy){
            mid <- mid + 1
            currSdy <- temp$study_accession[mid]
          }
          if(batch == 1){
            start <- 1
            end <- mid
          }else{
            start <- mid + 1
            end <- dim(temp)[1]
          }
          temp <- temp[start:end, ]
        }

        file_link <- paste0(
          self$config$labkey.url.base,
          "/_webdav/Studies/",
          temp$study_accession,
          "/%40files/rawdata/",
          folder,
          "/",
          URLencode(temp$file_info_name)
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
    }

    ..checkLinksOtherFolder <- function(folders, subdir, file_names){
      file_link <- paste0(
        self$config$labkey.url.base,
        "/_webdav/Studies/",
        folders,
        "/%40files/",
        subdir,
        file_names
      )

      file_exists <- unlist(
        mclapply(
          file_link,
          private$.checkUrl,
          mc.cores = mc.cores
        )
      )

      res <- data.frame(
        file_info_name = file_names,
        study_accession = folders,
        file_link = file_link,
        file_exists = file_exists,
        stringsAsFactors = FALSE
      )
    }
    # -----------------------------


    ## ------- MAIN ---------------
    startTime <- Sys.time()

    rawFolderData <- c("gene_expression_files", "fcs_sample_files", "fcs_control_files")

    if (file_type %in% rawFolderData) {
      folder <- ifelse(file_type == "gene_expression_files",
                       "gene_expression",
                       "flow_cytometry")
      res <- ..checkLinksRawFolder(file_type = file_type,
                                   folder = folder,
                                   batch = batch)

    } else {
      if (file_type == "protocols") {
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

        subdir <- "protocols/"
        file_names <- paste0(folders, "_protocol.zip")

      } else if (file_type == "gene_expression_matrices"){
        suppressWarnings(
          mx <- ImmuneSpaceR:::.getLKtbl(
            con = self,
            schema = "assay.ExpressionMatrix.matrix",
            query = "Runs",
            colNameOpt = "rname"
          )
        )

        folders <- mx$folder_name
        subdir <- "/analysis/exprs_matrices/"
        file_names <- paste0(mx$name, ".tsv")
      }

      res <- ..checkLinksOtherFolder(folders = folders,
                                     subdir = subdir,
                                     file_names = file_names)
    }

    ..messageResults(file_type = file_type,
                     file_exists = res$file_exists)

    endTime <- Sys.time()
    diff <- endTime - startTime
    message(diff, " ", attributes(diff)$units)

    res
  }
  # --------------------------
)

# Returns a string that can be used as a shell command on RServe machines
# (rsT / rsP) for downloading files for studies that are missing files
# according to checkRawFiles().
ISM$set(
  which = "public",
  name = "generateRawFilesCmd",
  value = function(rawFilesOutput) {
    sdys <- lapply(rawFilesOutput, function(x){
      x <- x[!x$file_exists,]
      return(unique(x$study_accession))
    })
    toDLstr <- paste(unique(unlist(sdys)), collapse = " ")
    if(toDLstr == ""){
      message("No missing files")
      return()
    }else{
      toDLstr <- paste("./getRawFiles.sh -sv ", toDLstr)
      return(toDLstr)
    }
  }
)


# PRIVATE ----------------------------------------------------------------------

# Check if the url exists (is accessible)
ISM$set(
  which = "private",
  name = ".checkUrl",
  value = function(url) {
    opts <- self$config$curlOptions

    res <- HEAD(url, config = opts)

    if (http_error(res)) {
      ret <- FALSE
    } else {
      if (http_type(res) == "application/json") {
        res <- httr::GET(url, config = opts)
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
