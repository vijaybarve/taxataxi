#' @importFrom bdvis gettaxo
#' @importFrom taxize gnr_resolve get_tsn get_colid synonyms synonyms_df get_tpsid

getIconicTab <- function(inatgr){
  iconic_taxon_name <- NULL
  Scientific_name <- NULL
  grpnames <- NULL
  inatgr <- as.data.frame(dplyr::summarise(dplyr::group_by(inatgr,iconic_taxon_name),
                                           ct=dplyr::n_distinct(Scientific_name)))
  inatgr1 <- inatgr
  inatgr$iconic_taxon_name <- grpnames$nps[match(inatgr$iconic_taxon_name, grpnames$inat)]
  i = which(is.na(inatgr[,1] == TRUE))
  inatgr$iconic_taxon_name <- as.character(inatgr$iconic_taxon_name)
  inatgr$iconic_taxon_name[i] = as.character(inatgr1$iconic_taxon_name[i])
  return(inatgr)
}

getClassasIconic <- function(inat){
  inat1 <- inat[which(is.na(inat$iconic_taxon_name)),]
  #print(paste("iconoc",dim(inat1)))
  nonmatch_taxo <- gettaxo(inat1,genus = TRUE,verbose = F,progress = F)
  if(dim(inat1)[1]>0){
    for(i in 1:dim(nonmatch_taxo)[1]){
      if(is.na(nonmatch_taxo$Class[i]) & !is.na(nonmatch_taxo$Phylum[i])){
        nonmatch_taxo$Class[i] <- nonmatch_taxo$Phylum[i]
      }
      inat$iconic_taxon_name[which(inat$Scientific_name==nonmatch_taxo$Scientific_name[i])] <- nonmatch_taxo$Class[i]

    }
  }
  return(inat)
}

getwords <- function(x,n) {
  return(paste(unlist(strsplit(x, split = "\\s+"))[1:n],collapse=" "))
}

toBinomial <- function(dat){
  for(i in 1:dim(dat)[1]){
    if(!is.na(dat$ScientificName[i])){
      dat$ScientificName[i] <- getwords(as.character(dat$ScientificName[i]),2)
    }
  }
  return(dat)
}

toBinomial_inat <- function(dat){
  for(i in 1:dim(dat)[1]){
    dat$Scientific_name[i] <- getwords(as.character(dat$Scientific_name[i]),2)
  }
  return(dat)
}


getInatRec <- function(namelist,inat_1p){
  dat <- inat_1p[inat_1p$Scientific_name %in% namelist,]
  return(dat)
}

list_itis_syn <- function(scname){
  #  syn <- synonyms(get_tsn_(scname, rows=1),db="itis")
  t1 <- NULL
  tsn <- get_tsn(scname, rows=1)[1]
  if(!is.na(tsn)){
    syn <- synonyms(tsn,db="itis")
    eval(parse(text=paste("t <- (syn$'",tsn,"')",sep='')))
    if(!is.null(t)){
      eval(parse(text=paste("t <- (syn$'",tsn,"'$acc_name)",sep='')))
      eval(parse(text=paste("t1 <- (syn$'",tsn,"'$syn_name)",sep='')))
      return(unique(t1))
    } else {
      return(NULL)
    }
  } else {
    return(NULL)
  }

}

#list_itis_syn(scname)
list_col_syn <- function(scname){
  #  syn <- synonyms(get_tsn_(scname, rows=1),db="itis")
  t1 <- NULL
  tsn <- get_colid(scname, rows=1)[1]
  if(!is.na(tsn)){
    syn <- synonyms(tsn,db="col")
    eval(parse(text=paste("t <- (syn$'",tsn,"')",sep='')))
    if(!is.null(t)){
      eval(parse(text=paste("t <- (syn$'",tsn,"'$acc_name)",sep='')))
      eval(parse(text=paste("t1 <- (syn$'",tsn,"'$name)",sep='')))
      return(unique(t1))
    } else {
      return(NULL)
    }
  } else {
    return(NULL)
  }
}

list_tropicos_syn <- function(scname){
  #  syn <- synonyms(get_tsn_(scname, rows=1),db="itis")
  t1 <- NULL
  tsn <- get_tpsid(scname, rows=1)[1]
  if(!is.na(tsn)){
    syn <- synonyms(tsn,db="tropicos")
    eval(parse(text=paste("t <- (syn$'",tsn,"')",sep='')))
    if(!is.null(t)){
      eval(parse(text=paste("t <- (syn$'",tsn,"'$acc_name)",sep='')))
      eval(parse(text=paste("t1 <- (syn$'",tsn,"'$scientificname)",sep='')))
      return(unique(t1))
    } else {
      return(NULL)
    }
  } else {
    return(NULL)
  }
}



get_itis_syn <- function(namelist){
  retset <- NULL
  for (i in 1:length(namelist)){
    set1 <- list_itis_syn(namelist[i])
    if(!is.null(set1)){
      set1 <- cbind(namelist[i],set1)
      retset <- rbind(retset,set1)
    }
  }
  return(retset)
}

get_all_syn <- function(namelist){
  retset <- NULL
  for (i in 1:length(namelist)){
    cat(paste(i,"of",length(namelist),"\n"))
    set1 <- list_itis_syn(namelist[i])
    set2 <- list_col_syn(namelist[i])
    set3 <- list_tropicos_syn(namelist[i])
    if(!is.null(set1)){
      set1 <- cbind(namelist[i],set1)
      retset <- rbind(retset,set1)
    }
    if(!is.null(set2)){
      set1 <- cbind(namelist[i],set2)
      retset <- rbind(retset,set1)
    }
    if(!is.null(set3)){
      if (set3 != "no syns found"){
        set1 <- cbind(namelist[i],set3)
        retset <- rbind(retset,set1)
      }
    }
  }
  #retset <- retset[which(retset$),]
  return(retset)
}

#get_itis_syn(namelist[1:10])
get_itis_binom_syn <- function(namelist){
  dat <- get_itis_syn(namelist)
  dat <- as.data.frame(dat)
  dat[] <- lapply(dat,as.character)
  names(dat) <- c("Acc_name","Syn")
  for(i in 1:dim(dat)[1]){
    dat$Syn[i] <- getwords(as.character(dat$Syn[i]),2)
  }
  return(unique(dat))
}

get_all_binom_syn <- function(namelist){
  dat <- get_all_syn(namelist)
  if(!is.null(dat)){
    dat <- as.data.frame(dat)
    dat[] <- lapply(dat,as.character)
    names(dat) <- c("Acc_name","Syn")
    for(i in 1:dim(dat)[1]){
      dat$Syn[i] <- getwords(as.character(dat$Syn[i]),2)
    }
    return(unique(dat))
  } else {
    return(NULL)
  }
}

#get_all_binom_syn(namelist[1:10])

get_binom_syn <- function(namelist){
  listdb <- c("itis","col")
  fulldat <- NULL
  dat2 <-
    for (i in 1:length(namelist)){
      dat1 <- c()
      for(j in 1:length(listdb)){
        dat <- synonyms(namelist[i],db=listdb[j])
        #dat1 <- synonyms_df(dat)
        dat1 <- c(dat1,dat)
      }
      dat2 <- synonyms_df(dat1)
      fulldat <- rbind(fulldat,dat1)
    }
  return((fulldat))
}

# Clean NPS Synonyms field
# Synonym field has some ids and formatting tags embeded

cleanSyn <- function(dat){
  dat$Synonyms <- gsub('[0-9]+', '', dat$Synonyms)
  #dat$Synonyms <- gsub('\\\\','',dat$Synonyms )
  dat$Synonyms <- gsub('<em>|</em>',' ',dat$Synonyms )
  #dat$Synonyms[which(dat$Synonyms !="")]
  dat$Synonyms <- gsub("(?<=[\\s])\\s*|^\\s+|\\s+$", "", dat$Synonyms, perl=TRUE)
  synbin <- NULL
  if(!all(is.na(dat$Synonyms))){
    for(i in 1:dim(dat)[1]){
      if(dat$Synonyms[i]!=""){
        t1 <- gsub('([[:upper:]])', '#\\1', dat$Synonyms[i])
        t2 <- trimws(strsplit(t1,"#")[[1]])
        if(t2[1]==""){
          t2 <- t2[-1]
        }
        t3 <- unique(t2)
        #print(t3)
        for(j in 1:length(t3)){
          synbin <- rbind(synbin,t3[j])
        }
        #print(t3)
      }
    }
    rownames(synbin) <- c()
    synbin <- as.data.frame(synbin)
    names(synbin) <- "ScientificName"
  }
  return(synbin)
}
