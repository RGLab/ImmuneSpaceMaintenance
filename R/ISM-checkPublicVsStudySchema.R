# Compare immport.public and study schema to determine studies that should have
# data copied over from ImmPort to Study via Immport Module "Copy Dataset"
# and data that has been changed in Immport that should be investigated
# prior to copying over.
#' @importFrom Rlabkey makeFilter
ISM$set(
  which = "public",
  name = "checkPublicVsStudySchema",
  value = function() {
    # Data from Immport schema (shown in study overview before login)
    pub <- labkey.selectRows(
      baseUrl = self$config$labkey.url.base,
      folderPath = "/Studies/",
      schemaName = "immport.public",
      queryName = "dimstudyassay",
      colSelect = c("Study", "Label"),
      colFilter = makeFilter(c(
        "Categorylabel",
        "EQUAL",
        "Raw data files"
      ))
    )

    # Data in study schema (shown in study overview after login)
    base <- "SELECT DISTINCT tmp.study_accession FROM tmp"

    fcsSampleFiles <- labkey.executeSql(
      baseUrl = self$config$labkey.url.base,
      folderPath = "/Studies/",
      schemaName = "study",
      sql = gsub("tmp", "fcs_sample_files", base)
    )
    fcsSampleFiles$Label <- "FCS sample files"

    fcsControlFiles <- labkey.executeSql(
      baseUrl = self$config$labkey.url.base,
      folderPath = "/Studies/",
      schemaName = "study",
      sql = gsub("tmp", "fcs_control_files", base)
    )
    fcsControlFiles$Label <- "FCS control files"

    geFiles <- labkey.executeSql(
      baseUrl = self$config$labkey.url.base,
      folderPath = "/Studies/",
      schemaName = "study",
      sql = gsub("tmp", "gene_expression_files", base)
    )
    geFiles$Label <- "Gene expression microarray data files"

    priv <- rbind(fcsControlFiles, fcsSampleFiles, geFiles)

    # Comparing data sets
    pub <- pub[pub$Study %in% unique(priv$`Study Accession`), ]
    pub$key <- paste(pub$Study, pub$Label)
    priv$key <- paste(priv$`Study Accession`, priv$Label)
    inPubNotPriv <- setdiff(pub$key, priv$key)
    inPrivNotPub <- setdiff(priv$key, pub$key)

    # return
    ret <- list(
      `In Public Not Study` = inPubNotPriv,
      `In Study Not Public` = inPrivNotPub
    )

    if (any(lengths(ret)) > 0) {
      message("Schema mismatches found. Check 'Publish to Study'.")
    }

    ret
  }
)
