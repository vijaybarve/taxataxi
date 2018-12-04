#' Build a iNaturalist data summary table with Park codes and counts of records
#' @param inat Clean iNaturalist data typically returned from \code{clean_inat_data}
#' @return data frame containing iNaturalist data summary with \code{unit_code}
#' and \code{count}
#' @importFrom stats complete.cases
#' @importFrom plyr count
#' @examples \dontrun{
#' inatsum <- inat_summary(inat)
#' }
#' @family data preparation
#' @export
inat_summary <- function(inat=NULL){
  if(is.null(inat)){
    stop("Please pass iNaturlaist data. Parks would be extracted form the data to build NPSecis master table")
  }
  outpath <- ""
  #parkcts <- sqldf("select unit_code,count(*) as ct from inat2016ip group by unit_code")
  parkcts <- count(inat,"Unit_code")
  parkcts <- parkcts[complete.cases(parkcts),]
  return(parkcts)
}
