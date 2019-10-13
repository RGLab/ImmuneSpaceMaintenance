#' @importFrom Rlabkey labkey.getFolders
#' @importFrom jsonlite fromJSON
#' @export
.listFiles <- function(link){
  response <- NULL
  res <- tryCatch(
    Rlabkey:::labkey.get(link),
    warning = function(w) return(w),
    error = function(e) return(NULL)
  )
  if (!is.null(res)) {
    tmp <- jsonlite::fromJSON(res, simplifyDataFrame = FALSE)
    response <- sapply(tmp$files, function(x) {
      return(x$text)
    }) # basename only
  }
  response
}
