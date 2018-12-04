#' Initial processing of iNaturalist project data
#'
#' Script to initiate data processing of iNaturalist project data
#' Uses downloaded .csv file of the iNaturalist project and nps_boundary shape file
#' Cleans the data to remove records from other years, missing geo-locations,
#' records outside parks
#' @param obsfile Full path to downloaded .csv file of the iNaturalist project
#' @param npsboundary filename with full path to NPS boundary shape file
#' @return data frame containing cleaned iNaturalist data
#' @importFrom utils write.csv
#' @importFrom readr read_csv
#' @examples \dontrun{
#' inat <- clean_inat_data(".\\data\\observations-17205.csv")
#' }
#' @family data preparation
#' @export
clean_inat_data <- function(obsfile,npsboundary=NA){
  # Read the downloaded file
  if (file.exists(obsfile)) {
    obs <- read_csv(obsfile)
  } else {
    stop("File",obsfile,"not found. Please specify filenems with full path")
  }

  # Format the data for bdvis
  inat <- bdvis::format_bdvis(obs,
                       Latitude="latitude",
                       Longitude="longitude",
                       Date_collected="observed_on",
                       Scientific_name="scientific_name")
  rm(obs)
  # Remove records form years other than 2016
  inat2016 <- inat[which(inat$Date_collected<as.Date("2017-01-01") & inat$Date_collected>as.Date("2015-12-31")),]
  rm(inat)
  # remove non geocoded records
  inat2016g <- inat2016[which(!is.na(inat2016$Latitude)),]
  sp::coordinates(inat2016g) <- ~ Longitude + Latitude
  rm(inat2016)

  # Read parks shape file
  if(is.na(npsboundary)){
    npsboundary <- ".\\data\\nps_boundary.shp"
  }
  parks <- maptools::readShapePoly(npsboundary)

  tmpinparks <- sp::over(inat2016g,parks)
  inat2016g <- as.data.frame(inat2016g)
  inat2016p <- cbind(inat2016g,tmpinparks$UNIT_CODE,tmpinparks$UNIT_NAME,tmpinparks$UNIT_TYPE)
  rm(inat2016g)

  # Rename matched fields
  if ("tmpinparks$UNIT_CODE" %in% names(inat2016p)) {
    names(inat2016p)[which(names(inat2016p)=="tmpinparks$UNIT_CODE")] <- "Unit_code"
    names(inat2016p)[which(names(inat2016p)=="tmpinparks$UNIT_NAME")] <- "Unit_name"
    names(inat2016p)[which(names(inat2016p)=="tmpinparks$UNIT_TYPE")] <- "Unit_type"
  }

  inat2016ip <- inat2016p[which(!is.na(inat2016p$Unit_code)),]
  rm(inat2016p)
  write.csv(inat2016ip,"inat2016ip.csv",row.names = F )
  return(inat2016ip)
}
