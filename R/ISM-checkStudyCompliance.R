#' @include ISM.R
NULL



# PUBLIC -----------------------------------------------------------------------

# This method allows admin to check which studies are compliant with the
# following modules/data files (GEF, RAW, GEO, GEM, DE, GEE, IRP, GSEA, DR)
# The return value can be either a dataframe or a summary list (format option).
# The summaryList is a list of only non-compliant studies, with summary information.
#' @importFrom Rlabkey labkey.selectRows labkey.executeSql
ISM$set(
  which = "public",
  name = "checkStudyCompliance",
  value = function(reload = FALSE,
                   summarize = TRUE,
                   filterNonGE = TRUE,
                   showAllCols = FALSE,
                   onlyShowNonCompliant = TRUE,
                   verbose = FALSE) {
    ## ----- HELPERS ----------
    # Get list of studies turned on for a specific module
    ..getModSdys <- function(name) {
      url <- url <- paste0(baseUrl, "/immport/studies/containersformodule.api?name=", name)
      res <- unlist(lapply(rjson::fromJSON(Rlabkey:::labkey.get(url))[[1]], function(x) {
        x[["name"]]
      }))
      res <- .spSort(res[grepl("SDY[0-9]+", res)])
    }

    ## ------- MAIN -----------
    # Check cache first and use preferentially
    if (!is.null(self$cache[["complianceDF"]]) && !reload) {
      compDF <- self$cache[["complianceDF"]]
    } else {
      baseUrl <- self$config$labkey.url.base # For labkey.executeSql calls

      # Prepare list of studies for generating results table rownames
      if (self$study != "Studies") {
        Sdys <- c(self$study) # single study
      } else {
        Sdys <- private$.getSdyVec()
        Sdys <- .spSort(Sdys)
      }

      # Prepare list of columns included in the table from Modules and Gene-Expression work
      mods <- c(
        # ---- Gene Expression Data ----
        "GEF", # Gene Expression Files query for meta-data
        "RAW", # Raw data on RServe filesystem and imported from ImmPort
        "GEO", # Raw data in Gene Expression Omnibus Dbase
        "GEM_implied", # Gene Expression Matrix
        "GEM_actual",
        # ---- Modules ----
        "DE_implied", # Data Explorer
        "DE_actual",
        "GEE_implied", # Gene Expression Explorer
        "GEE_actual",
        "DGEA_implied", # Differential Expression Analysis (report)
        "DGEA_actual",
        "DGEA_missing",
        "IRP_implied", # Immune Response Predictor
        "IRP_actual",
        "IrpTimepoints",
        "GSEA_implied", # Gene Set Enrichment Analysis
        "GSEA_actual",
        "DR_implied", # Dimension Reduction
        "DR_actual"
      )

      compDF <- data.frame(
        matrix(
          nrow = length(Sdys),
          ncol = length(mods)
        ),
        row.names = Sdys
      )
      colnames(compDF) <- mods

      # Get lists of studies currently enabled for each module
      compDF$DE_actual <- rownames(compDF) %in% ..getModSdys("DataExplorer")
      compDF$GEE_actual <- rownames(compDF) %in% ..getModSdys("GeneExpressionExplorer")
      compDF$GSEA_actual <- rownames(compDF) %in% ..getModSdys("GeneSetEnrichmentAnalysis")
      compDF$IRP_actual <- rownames(compDF) %in% ..getModSdys("ImmuneResponsePredictor")
      compDF$DGEA_actual <- rownames(compDF) %in% ..getModSdys("DifferentialExpressionAnalysis")
      compDF$DR_actual <- rownames(compDF) %in% ..getModSdys("DimensionReduction")

      ################################
      ###   Gene Expression Data   ###
      ################################

      # GEF
      # --------
      # study.gene_expression_files: table with metadata about gene expresssion data from immport
      gef <- self$getDataset("gene_expression_files")
      compDF$GEF <- rownames(compDF) %in% .subidsToSdy(gef$participant_id)

      # RAW
      # --------
      # information on how to find the raw gene expression data flat files on Rserve machine
      file_list <- private$.getGEFileNames(TRUE)
      file_list <- file_list[file_list != "NULL"]
      compDF$RAW <- rownames(compDF) %in% names(file_list)[file_list == TRUE]

      # GEO
      # --------
      # geo accession IDs for raw data in gene expression omnibus database
      geoGef <- gef[!is.na(gef$geo_accession), ]
      compDF$GEO <- rownames(compDF) %in% .subidsToSdy(geoGef$participant_id)

      # GEM_implied
      # --------
      # Gene Expression Matrix can be created from raw files on Rserve or files downloaded from GEO.
      compDF$GEM_implied <- compDF$RAW == T | compDF$GEO == T

      # GEM_actual
      # --------
      # Does study have gene expression matrices in the assay.ExpressionMatrix.matrix.Runs query?
      studiesWithGems <- unique(self$cache$GE_matrices$folder)
      compDF$GEM_actual <- rownames(compDF) %in% studiesWithGems


      ################################
      ###          Modules         ###
      ################################

      # Get immune response data for GEE and IRP
      # `hai` is hemaglutinin Inhibition assay data and `nab` is neutralizing antibody titer
      # assay data.  Both assays measure the antibodies present for specific viral antigens.
      # HAI is predominantly used for Influenza studies while NAb is more common for others.
      hai <- self$getDataset("hai")
      nab <- self$getDataset("neut_ab_titer")
      immuneResponse <- rbind(nab, hai, fill = TRUE) # nab has a col that hai does not


      # GEE - Gene Expression Explorer
      # --------
      # visualization of expression level of 1 or more genes
      # vs immune response (HAI or NAb):
      # Must have subjects with both GEM from any timepoint and response data for any timepoint
      # (timepoints do not have to be the same)
      # NOTE: when GEE is changed to allow NAb, can uncomment respSubs lines
      # resp <- union(resp$participant_id, nab$participant_id)
      inputSmpls <- labkey.selectRows(
        baseUrl = baseUrl,
        folderPath = "/Studies",
        schemaName = "study",
        queryName = "HM_InputSamplesQuery",
        containerFilter = "CurrentAndSubfolders",
        colNameOpt = "rname"
      )
      setDT(inputSmpls)

      inputSmpls$study <- gsub("SUB[^>]+\\.", "SDY", inputSmpls$participantid)

      exprResp <- merge(
        inputSmpls,
        hai,
        by.x = c("participantid", "study_time_collected"),
        by.y = c("participant_id", "study_time_collected")
      )

      compDF$GEE_implied <- rownames(compDF) %in% .subidsToSdy(unique(exprResp$participantid))

      # IRP - Immune Response Predictor
      # --------
      # This module can be used to automatically select a group
      # of genes whose expression at a given time point (e.g. gene expression levels at day 0)
      # best predicts a given immunological response at a later time point (e.g. HAI at day 28)
      # Requires studies with subjects from multiple cohorts with GEM data at both target
      # timepoint and baseline + response data that has baseline and a later timepoint.
      # GEM and response later timepoints do NOT need to be the same!

      resp <- immuneResponse[, list(study_time_collected,
                                    study_time_collected_unit,
                                    response = value_preferred / mean(value_preferred[study_time_collected <= 0], na.rm = TRUE)
                                    ),
                                by = "virus,participant_id"
                             ]
      resp <- resp[ !is.na(response) ]

      # NOTE: At least SDY180 has overlapping study_time_collected for both hours and days
      # so it is important to group by study_time_collected_unit as well. This is reflected
      # in IRP_timepoints_hai/nab.sql.

      # Subset to only participants with response data
      geCohortSubs <- inputSmpls[ participantid %in% resp$participant_id ]
      # Subset to only samples from studies where there is data from multiple cohorts at
      # a given timepoint
      geCohortSubs <- geCohortSubs[, .SD[length(unique(cohort)) > 1],
                                     by = .(study, study_time_collected, study_time_collected_unit)
                                   ]
      # Subset to only samples where there is baseline data and data from other timepoints
      geCohortSubs <- geCohortSubs[, .SD[length(unique(study_time_collected)) > 1 & 0 %in% unique(study_time_collected)],
                                     by = .(study, cohort, study_time_collected_unit)
                                   ]
      compDF$IRP_implied <- rownames(compDF) %in% unique(geCohortSubs$study)

      # Get IrpTimepoints
      # TODO:  Change to "IRP_missing" to be consistent, and only include when noncompliant (or missing?)
      # TODO:  Determine if this field is necessary
      studyTimepoints <- geCohortSubs[, list(timepoints = paste(sort(unique(study_time_collected)),
                                                                collapse = ",")),
                                        by = .(study)
                                      ]
      compDF$IrpTimepoints <- studyTimepoints$timepoints[ match(rownames(compDF), studyTimepoints$study) ]

      # DGEA - Differential Expression Analysis
      # --------
      # creates GEAR and GEA tables
      # compares baseline and other timepoint expression levels to find
      # genes that are sigificantly differentially expressed between baseline
      # and other timepoints

      # Has DGEA already been run? If so, if data has been added gea might not be complete
      # with all timepoints and needs to be rerun
      # Check: Do gea results exist? are they complete?
      existGEA <- labkey.selectRows(
        baseUrl = baseUrl,
        folderPath = "/Studies/",
        schemaName = "gene_expression",
        queryName = "gene_expression_analysis",
        colNameOpt = "rname",
        showHidden = TRUE
      )

      containers <- labkey.selectRows(
        baseUrl = baseUrl,
        folderPath = "/Studies/",
        schemaName = "core",
        queryName = "Containers",
        containerFilter = "CurrentAndSubfolders",
        showHidden = TRUE
      )

      existGEA$sdy <- containers$`Display Name`[ match(existGEA$container, containers$`Entity Id`)]

      gea <- suppressWarnings(lapply(studiesWithGems, FUN = function(sdy) {
        # TODO:  use inputSmpls instead and subset
        impliedGEA <- inputSmpls[study == sdy]

        # ---- summarize to have same form and info as currGEA ----

        # 1. Remove all arm_name * study_time_collected with less than 4 replicates
        # otherwise predictive modeling cannot work
        impliedGEA[, subs := length(unique(participantid)),
                     by = .(cohort, study_time_collected, study_time_collected_unit)
                   ]
        impliedGEA <- impliedGEA[ subs > 3 ]

        # 2. Check for baseline within each arm_name and then filter out baseline
        impliedGEA[, baseline := any(study_time_collected <= 0), by = .(cohort) ]
        impliedGEA <- impliedGEA[ baseline == TRUE ] # filter out arms with no baseline
        impliedGEA <- impliedGEA[ study_time_collected > 0 ] # remove baseline

        # 3. Generate key

        impliedGEA[, key := paste(cohort_type, study_time_collected, study_time_collected_unit)]

        # 4. Summarize by arm_name * study_time_collected for number of subs and key
        smryGEA <- impliedGEA[, list(key = unique(key), subs = unique(subs)),
                                by = .(cohort_type, study_time_collected, study_time_collected_unit)
                              ]

        # -------------------------------------------

        # Get current GEA and compare
        currGEA <- existGEA[ existGEA$sdy == sdy, ]
        currGEA$key <- paste(currGEA$arm_name, currGEA$coefficient)

        if (nrow(smryGEA) > 0) {
          diff <- sort(setdiff(smryGEA$key, currGEA$key)) # In implied and NOT in current
          missing_data <- if (length(diff) == 0) {
            "no diff"
          } else {
            paste(diff, collapse = "; ")
          }
        } else {
          missing_data <- NA
        }

        res <- c("DGEA_implied" = nrow(smryGEA) > 0, "DGEA_missing" = missing_data)
      }))

      names(gea) <- studiesWithGems

      gea <- data.frame(do.call(rbind, gea), stringsAsFactors = FALSE)
      compDF[studiesWithGems, "DGEA_implied"] <- gea$DGEA_implied
      compDF[studiesWithGems, "DGEA_missing"] <- gea$DGEA_missing
      compDF$DGEA_implied[ is.na(compDF$DGEA_implied) ] <- FALSE
      compDF$DGEA_implied <- as.logical(compDF$DGEA_implied)

      # GSEA - Gene set enrichment analysis
      # --------
      # visualize how gene expression changes over time
      # for groups of genes:
      # studies with subjects having results in the
      # gene_expression.gene_expression_analysis_results (GEAR) table for multiple non-baseline timepoints
      # (GEA and GEAR are generated by differential expression analysis (DGEA) report)
      gearSql <- "SELECT DISTINCT analysis_accession.coefficient FROM gene_expression_analysis_results"
      gear <- sapply(studiesWithGems, FUN = function(sdy) {
        res <- suppressWarnings(
          tryCatch(
            labkey.executeSql(
              baseUrl = baseUrl,
              folderPath = paste0("/Studies/", sdy),
              schemaName = "gene_expression",
              sql = gearSql
            ),
            error = function(e) {
              return(NA)
            }
          )
        )
        output <- !is.na(res) && nrow(res) > 1
      })

      compDF$GSEA_implied <- rownames(compDF) %in% names(gear)[gear == TRUE]

      # DE - Data Explorer
      # --------
      # visualize assay data b/c ISC_study_datasets
      # cannot provide gene_expression info, we use compDF$DGEA_actual as a
      # proxy since it pulls the current GEA query, which should have the same
      # info as DGEA_filteredGEAR ( what the con$plot() uses via
      # con$getGEAnalysis() )
      deSets <- c(
        "Neutralizing antibody titer",
        "Enzyme-linked immunosorbent assay (ELISA)",
        "Enzyme-Linked ImmunoSpot (ELISPOT)",
        "Hemagglutination inhibition (HAI)",
        "Polymerisation chain reaction (PCR)",
        "Flow cytometry analyzed results",
        "Multiplex bead array asssay"
      )
      compDF$DE_implied <- sapply(rownames(compDF), FUN = function(sdy) {
        res <- suppressWarnings(
          tryCatch(
            labkey.executeSql(
              baseUrl = baseUrl,
              folderPath = paste0("/Studies/", sdy),
              schemaName = "study",
              sql = "SELECT Label FROM ISC_datasets"
            ),
            error = function(e) {
              return(NA)
            }
          )
        )
        ret <- any(res[[1]] %in% deSets) | compDF$DGEA_actual[rownames(compDF) == sdy]
      })

      # DR - Dimension Reduction (DR)
      # --------
      # visualize multiple assays/features using dimension reduction algorithms
      # to identify clustering. Requires at enough subjects and features to come
      # up with at least a 3x3 matrix. See note below for details on how this is
      # determined.
      #
      # NOTE:  For now, it seems that simply checking to see if there are any
      # assay/subject combinations where number of subjects and features are
      # both greater than 3 is a good enough estimate of whether or not
      # dimension reduction makes sense or is possible without doing too much of
      # the checks and filtering that already happens in the module. It may mark
      # as false some studies where dimension reduction might be possible when
      # including multiple timepoints or assays. If this proves to be a problem,
      # we could further group by timepoint or assay to get an idea of what the
      # dimensions would be.

      # get dimension reduction assay info
      dimRedux_assay_data <- labkey.selectRows(
        baseUrl = baseUrl,
        folderPath = "/Studies",
        schemaName = "study",
        queryName = "DimRedux_assay_data_computed",
        containerFilter = "CurrentAndSubfolders",
        colNameOpt = "rname"
      )

      # Add a column for study
      dimRedux_assay_data$study <- gsub("SUB[^>]+\\.", "SDY", dimRedux_assay_data$participantid)
      setDT(dimRedux_assay_data)

      # Group by study, timepoint, and assay, and get the number of subjects and features for that
      # assay
      dimensionInfo <- dimRedux_assay_data[, .(
        subjectCount = length(unique(participantid)),
        featureCount = min(features)
      ),
      by = c("study", "timepoint", "name")
      ]

      # Are there any lines where subject count and feature count are both greater than three?
      dimensionInfo[, dimMinMet := subjectCount >= 3 & featureCount >= 3]
      dimRedPossible <- dimensionInfo[, .(dimMinMet = any(dimMinMet)), by = "study"]

      compDF$DR_implied <- rownames(compDF) %in% dimRedPossible[dimMinMet == TRUE, study]



      colOrder <- c(
        "RAW",
        "GEF",
        "GEO",
        "GEM_implied",
        "GEM_actual",
        "DE_implied",
        "DE_actual",
        "GEE_implied",
        "GEE_actual",
        "DGEA_implied",
        "DGEA_actual",
        "DGEA_missing",
        "IRP_implied",
        "IRP_actual",
        "IrpTimepoints",
        "GSEA_implied",
        "GSEA_actual",
        "DR_implied",
        "DR_actual"
      )

      rowOrder <- Sdys[order(
        gsub("([A-Z]+)([0-9]+)", "\\1", Sdys),
        as.numeric(gsub("([A-Z]+)([0-9]+)", "\\2", Sdys))
      )]

      compDF <- compDF[order(match(row.names(compDF), rowOrder)), order(match(colnames(compDF), colOrder))]

      # Cache ---------------
      self$cache[["complianceDF"]] <- compDF
    }

    ################################
    #  Studies that are loaded but not enabled
    ################################

    # Get list of study folders from webdav
    folder_link <- paste0(
      self$config$labkey.url.base,
      "/_webdav/Studies?method=JSON"
    )
    shareStudies <- grep("SDY\\d+", private$.listISFiles(folder_link), value = TRUE)

    # Get list of studies on IS
    conStudies <- labkey.selectRows(
      baseUrl = self$config$labkey.url.base,
      folderPath = "/home",
      schemaName = "lists",
      queryName = "Studies",
      viewName = "",
      colSort = "id",
      colFilter = NULL,
      containerFilter = NULL
    )

    missingStudies <- setdiff(shareStudies, conStudies$Name)
    if (length(missingStudies > 0)) {
      message(
        paste0(
          length(missingStudies), " studies present on webdav but not enabled: ",
          paste(missingStudies, collapse = ", ")
        )
      )
    }


    ################################
    ###     Filter/Summarize     ###
    ################################

    if (summarize) {

      # Get noncompliant studies
      modules <- c("GEM", "DE", "GEE", "IRP", "GSEA", "DGEA", "DR")
      compliant <- data.frame(lapply(modules, function(module) {
        impl <- grep(paste0(module, "_implied"), colnames(compDF))
        act <- grep(paste0(module, "_actual"), colnames(compDF))

        if (module == "DGEA") {
          imp_vs_act <- compDF[[impl]] == compDF[[act]]
          missing_dat <- is.na(compDF["DGEA_missing"]) | compDF["DGEA_missing"] == "no diff"
          return(compliant <- imp_vs_act == missing_dat)
        } else {
          return(compliant <- compDF[[impl]] == compDF[[act]])
        }
      }))

      colnames(compliant) <- modules
      rownames(compliant) <- rownames(compDF)
      nonCompliantStudies <- rownames(compliant[ !apply(compliant, 1, all), ])

      summaryList <- lapply(nonCompliantStudies, function(study) {
        sl <- list(
          modules = modules[!compliant[study, ]]
        )
        if ("IRP" %in% sl$modules) {
          sl$IrpTimepoints <- compDF[study, "IrpTimepoints"]
        }
        if ("DGEA" %in% sl$modules) {
          sl$DGEA_missing <- compDF[study, "DGEA_missing"]
        }
        return(sl)
      })
      names(summaryList) <- nonCompliantStudies
      return(summaryList)
    } else {

      # Filter out studies that don't have GE since this is basis for everything
      if (filterNonGE) {
        compDF <- compDF[compDF$GEO | compDF$GEF, ]
      }
      # Subset to only show problematic studies
      if (onlyShowNonCompliant) {
        redux <- compDF[, grep("implied|actual|missing", colnames(compDF))]
        mod_sub <- c("DE", "GEE", "IRP", "GSEA", "DGEA", "DR")
        compliant <- lapply(mod_sub, FUN = function(mod) {
          idx <- grepl(mod, names(redux))
          sub <- redux[, idx]
          if (mod == "DGEA") {
            imp_vs_act <- sub[, 1] == sub[, 2]
            missing_dat <- is.na(sub[, 3]) | sub[, 3] == "no diff"
            compliant <- imp_vs_act == missing_dat
          } else {
            compliant <- sub[, 1] == sub[, 2]
          }
        })
        compliant[[7]] <- compDF$GEM_implied == compDF$GEM_actual
        compliant <- do.call(cbind, compliant)
        row.names(compliant) <- row.names(redux)
        idx <- which(apply(compliant, 1, all))
        compDF <- compDF[-(idx), ]
      }

      # Defaults to showing only the actual module status and the difference with the implied
      if (!showAllCols) {
        compDF <- compDF[, grep("act|implied$", colnames(compDF))]
      }

      if (verbose) {
        message("NOTE: \n Return objects have an actual column that was generated by a call to the module url and an implied column that \n was created by looking at the filesystem.")
      }

      return(compDF)
    }
  }
)



# PRIVATE ----------------------------------------------------------------------



# HELPER -----------------------------------------------------------------------

# Sort studies by number
.spSort <- function(vec) {
  if (length(vec) > 0) {
    vec <- sort(as.numeric(gsub("SDY", "", vec)))
    vec <- paste0("SDY", as.character(vec))
  } else {
    vec <- NULL
  }

  vec
}


# Get SDY IDs from subids
.subidsToSdy <- function(subids) {
  sdys <- unique(gsub("^SUB.+", NA, unlist(strsplit(subids, split = "\\."))))
  sdys <- sdys[!is.na(sdys)]
  sdys <- paste0("SDY", sdys)

  sdys
}

