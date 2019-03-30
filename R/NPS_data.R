#' Download all individual Park checklist files
#'
#' Script to download individual Park checklist files
#' @param inatsum iNaturalist data summary typically returned from \code{inat_summary}
#' @param outpath Output path to store the NPSpecies files and results
#' @return data frame containing iNaturalist data summary with \code{unit_code}
#' and \code{count}
#' @importFrom utils read.csv
#' @examples \dontrun{
#' nps <- NPS_data(inatsum)
#' }
#' @family data preparation
#' @export
NPS_data <- function(inatsum=NULL,outpath=NA){
  if(is.null(inatsum)){
    stop("Please pass iNaturlaist data. Parks would be extracted frmm the data to build NPSecis master table")
  }
  if(is.na(outpath)){
    outpath <- ""
  }
  synfile <- ".\\data\\synlist.csv"
  if (file.exists(synfile)){
    synlist <- read.csv(synfile)
  } else {
    synlist <- NULL
  }
  parkcts <- inatsum
  for (i in 1:dim(parkcts)[1]){
    NPSd <- get_NPS_list(parkcts$Unit_code[i], "detaillist")
    filename <- paste(outpath,parkcts$Unit_code[i],"_det.csv",sep = '')
    if(dim(NPSd$dat)[1]>0){
      write.csv(NPSd$dat,filename,row.names = F)
    }
    ### If there are synonyms in the checklist add them to synlist
    if(!is.null(NPSd$Syn)){
      parksyn <- clean_NPS_syn(NPSd$Syn)
      print(dim(parksyn)[1])
      parksyn <- parksyn[which(!(parksyn$Acc_name %in% synlist$Acc_name)),]
      print(dim(parksyn)[1])
      synlist <- rbind(synlist,parksyn)
      write.csv(synlist,synfile,row.names = F)
      print(paste("Synonyms file stored:",synfile,"with records:",dim(parksyn)[1]))
    }
  }

  # Now merge these files to get a NPS data in single file

  retmat <- NULL
  for (i in 1:dim(parkcts)[1]){
    filename <- paste(outpath,parkcts$Unit_code[i],"_det.csv",sep = '')
    if(file.exists(filename)){
      spdat <- read.csv(filename)
      if(dim(spdat)[1]>0){
        print(i)
        #spdat1 <- spdat[,fldnames]
        retmat <- rbind.match.columns(retmat,spdat)
      }
    }
  }

  NPSall <- retmat
  outfile <- paste(outpath,"NPSfull.csv",sep='')
  write.csv (NPSall, outfile, row.names = F)
  print(paste("Output is stored in:",outfile))
  return(NPSall)
}

## Match columns and then return the matching data
rbind.match.columns <- function(input1, input2) {
  if(is.null(input2)){
    return(input1)
  }
  if(is.null(input1)){
    return(input2)
  }

  n.input1 <- ncol(input1)
  n.input2 <- ncol(input2)
  if (n.input2 < n.input1) {
    TF.names <- which(names(input2) %in% names(input1))
    column.names <- names(input2[, TF.names])
  } else {
    TF.names <- which(names(input1) %in% names(input2))
    column.names <- names(input1[, TF.names])
  }

  return(rbind(input1[, column.names], input2[, column.names]))
}

