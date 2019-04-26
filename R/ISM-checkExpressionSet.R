#' @include ISM.R
NULL



# PUBLIC -----------------------------------------------------------------------

# getGEMatrix test function
ISM$set(
  which = "public",
  name = "checkExpressionSet",
  value = function(allMatrices = FALSE,
                   ...) {
    # Grab expression set
    mat_names <- ifelse(allMatrices == FALSE, self$cache$GE_matrices$name[[1]], self$cache$GE_matrices$name)
    es <- self$getGEMatrix(mat_names, ...)
    opts <- list(...)

    if (length(opts) == 0) {
      opts <- list()
      opts$outputType <- "summary"
    } else {
      opts <- list(...)
    }
    # expression matrix +
    em <- Biobase::exprs(es)
    pd <- Biobase::pData(es)

    res <- cbind(.checkEM(em, opts, self), .checkPD(pd, self), .checkBiosample(em, pd))
    return(res)
  }
)



# PRIVATE ----------------------------------------------------------------------



# HELPER -----------------------------------------------------------------------

# Check expression matrix
.checkEM <- function(em, opts, self) {
  # get feature set
  anno <- labkey.selectRows(
    baseUrl = self$config$labkey.url.base,
    folderPath = self$config$labkey.url.path,
    schemaName = "Microarray",
    queryName = "FeatureAnnotation",
    colSelect = c("FeatureId", "GeneSymbol"),
    maxRows = 20,
    showHidden = TRUE
  )

  # compare to gem rows
  anno <- ifelse(opts$outputType == "summary", anno$`Gene Symbol`, anno$`Feature Id`)
  anno_match <- all(anno %in% row.names(em))
  # check range (log2)
  expr_within_range <- all(0 < range(em) & range(em) < 30)
  # check num of genes
  min_genes <- ifelse(opts$outputType == "summary", 10000, 20000)
  gene_num <- length(row.names(em)) >= min_genes
  res <- data.frame(outputType = opts$outputType, anno_match, expr_within_range, gene_num)
  return(res)
}


# Check pdata
.checkPD <- function(pd, self) {
  cohort_type_col <- "cohort_type" %in% colnames(pd)
  ct_split <- do.call(rbind, strsplit(pd$cohort_type, "_", fixed = TRUE))
  # does cohort type cohort match pd$cohort
  cohort_match <- all(ct_split[, 1] == pd$cohort)
  # does type == labkey lookup for cell type
  lk_smpl_type <- labkey.selectRows(
    baseUrl = self$config$labkey.url.base,
    folderPath = self$config$labkey.url.path,
    schemaName = "immport",
    queryName = "lk_sample_type",
    showHidden = TRUE,
    colSelect = "Name"
  )
  type_match <- all(ct_split[, 2] %in% lk_smpl_type$Name)
  res <- data.frame(cohort_type_col, cohort_match, type_match)
  return(res)
}


# Check biosamples from pdata and expression matrix
.checkBiosample <- function(em, pd) {
  # change to all.equal fxn
  biosample_match <- all.equal(row.names(pd), colnames(em))
  if (all(biosample_match != TRUE)) {
    biosample_match <- FALSE
  }
  res <- data.frame(biosample_match)
  return(res)
}
