#' Clean NPSpecies synonyms list
#'
#' @param npssyn Raw list of NPSpecies synonyms returned by \code{get_NPS_list}
#' @return list of clean synonyms in a dat frame
#' @examples \dontrun{
#' clean_NPS_syn(npssyn)
#' }
#' @family data preparation
#' @export
clean_NPS_syn <- function(npssyn=NULL){
  retdat <- NULL
  if(!is.null(npssyn)){
    npssyn$Synonym <- gsub('var.|spp.',' ',npssyn$Synonym )
    npssyn$Synonym <- gsub("^ *|(?<= ) | *$", "", npssyn$Synonym, perl = TRUE)
    retdat <- npssyn[,c("ScientificName","Synonym")]
    names(retdat) <- c("Acc_name","Syn")
  }
  return(retdat)
}
