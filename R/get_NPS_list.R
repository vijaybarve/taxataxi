#' Get NPSpecies data for single park
#'
#' @param unitcode Four letter Unit code for the park
#' @param listtype Type of list to be downloaded. Accepted values are "fulllist" for
#'   for complete list, "checklist" for checklist of species, "detaillist" for
#'   full list of species with details
#' @return list of two data frames first \code{data} containing list of species and
#'  second \code{Syn} containing list of synonyms with respective accepted names
#' @examples \dontrun{
#' getNPSlist(unitcode="ABLI", listtype="fulllist")
#' }
#' @family data preparation
#' @importFrom jsonlite fromJSON
#' @export
get_NPS_list <- function(unitcode = "ABLI", listtype="fulllist"){
  url <- paste("https://irmaservices.nps.gov/v3/rest/npspecies/",listtype,"/",unitcode,"/",sep="")
  print(url)
  doc <- as.data.frame(fromJSON(url))
  ### Process synonyms
  Syndat <- NULL
  for( i in 1:dim(doc)[1]){
    if(length(doc$Synonyms[i][[1]])>0){
      if (dim(doc$Synonyms[i][[1]])[1]>0){
        recs <- cbind(doc$Id[i],doc$TaxaCode[i], doc$ScientificName[i],doc$Synonyms[i][[1]])
        Syndat <- rbind(Syndat,recs)
      }
    }
  }
  if(!is.null(Syndat)){
    names(Syndat) <- c("Id","TaxaCode","ScientificName","SynTaxaCode","Synonym","SynFormated")
  }
  ### End Process Synonyms
  retdat <- NULL
  retdat$data <- doc[ , -which(names(doc) %in% c("Synonyms"))]
  retdat$Syn <- Syndat
  return(retdat)
}
