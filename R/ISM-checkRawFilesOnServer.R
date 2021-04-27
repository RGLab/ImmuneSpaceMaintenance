#' @include ISM.R
NULL

# PUBLIC -----------------------------------------------------------------------

ISM$set(
  which = "public",
  name = "checkRawFilesOnServer",
  value = function(file_type, summarizeByStudy = TRUE) {
    temp <- self$getDataset(file_type, original_view = TRUE)

    # Create local address
    if (file_type == "fcs_control_files") {
      temp <- temp[, file_info_name := control_file]
      temp <- temp[, c("pid", "sid") := data.table::tstrsplit(participant_id, "\\.")]
      temp <- temp[, study_accession := paste0("SDY", sid)]
    }

    temp <- temp[!is.na(file_info_name)]
    temp <- unique(temp[, list(study_accession, file_info_name)])
    temp$path <- paste0(
      "/share/files/Studies/",
      temp$study_accession,
      "/@files/rawdata/flow_cytometry/",
      temp$file_info_name
    )

    # Check in filesystem
    temp$file_exists <- unlist(parallel::mclapply(temp$path, file.exists, mc.cores = 4))
    temp <- temp[file_exists != TRUE, ]

    if (summarizeByStudy) {
      temp <- temp[, .(missing_files = .N), by = study_accession]
    }

    return(temp)
  }
)
