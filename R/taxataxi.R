#' main function
#'
#' Script to
#' @param unitcode unitcode of the park
#' @param rg If the function should consider only research grade
#' records (default FALSE)
#' @param inat iNaturalist data
#' @param NPSall NPS data
#' @return list
#' @examples \dontrun{
#' taxataxi("ABLI")
#' }
#' @family data preparation
#' @importFrom stats na.omit
#' @importFrom plyr rbind.fill
#' @importFrom sqldf sqldf
#' @export
taxataxi <- function(inat,NPSall,unitcode,rg=FALSE){
  parkstat <- NULL
  cat(paste("\n \nPark code",unitcode,"\n"))

  #GEt data for the park from iNaturalist project download and NPS full file
  inat_1p <- sqldf(paste("Select * from inat where unit_code = '",unitcode,"'",sep = ""))
  nps_1p <- sqldf(paste("Select * from NPSall where unitcode  = '",unitcode,"'",sep = ""))
  nps_1p[] <- lapply(nps_1p, as.character)
  if(rg){
    inat_1p <- inat_1p[which(inat_1p$quality_grade=="research"),]
  }

  # Creating new records in NPS for all the synonyms so matching is simple with iNat
  #tmpsyn <- cleanSyn(nps_1p)
  tmpsyn <- clean_NPS_syn(nps_1p)
  nps_1p_orig <- nps_1p
  nps_1p <- rbind.fill(nps_1p,tmpsyn)

  # Remove all names except binomials
  nps_1p <- toBinomial(nps_1p)

  ScientificName <- NULL
  Category <- NULL
  cat(paste("iNaturalist project has",dim(inat_1p)[1],"records.\n"))
  cat(paste("NPS has",dim(nps_1p)[1],"taxa.\n"))
  nps_1p_tab <- as.data.frame(dplyr::summarise(dplyr::group_by(nps_1p,Category),
                                                ct=dplyr::n_distinct(ScientificName)))
  print(nps_1p_tab)

  parkstat$parkcode <- unitcode
  parkstat$iNatRec <- dim(inat_1p)[1]
  parkstat$NPSrec <- dim(nps_1p)[1]
  parkstat$NPSrec_sum <- nps_1p_tab

  #Distinct taxa in iNaturalist
  # Need to figure out how can all records pulled out for subsets later
  inat_1p_taxa <- sqldf("select distinct taxon_species_name, iconic_taxon_name from inat_1p")
  inat_1p_taxa <- na.omit(inat_1p_taxa)
  #inat_1p_taxa <- toBinomial_inat(inat_1p_taxa)

  cat(paste("iNaturalist project has",dim(inat_1p_taxa)[1],"unique taxa.\n"))
  parkstat$inattaxa <- dim(inat_1p_taxa)[1]

  # Resolve iNaturalist names with ITIS
  inat_1p_itis <- as.data.frame(gnr_resolve(names=inat_1p_taxa$taxon_species_name, data_source_ids=3))
  # Get ITIS unresolved names
  inat_1p_noitis <- sqldf("select * from inat_1p_taxa where taxon_species_name not in (select user_supplied_name from inat_1p_itis)")
  inat_1p$itis <- TRUE
  inat_1p$itis[which(inat_1p$taxon_species_name %in% inat_1p_noitis$taxon_species_name)] <- FALSE

  # Resolve them with other taxonomies
  #inatname01 <- as.data.frame(gnr_resolve(names=inatnoitis$taxon_species_name))
  no_inat_ITIS_match <- sqldf("select count(distinct user_supplied_name) from inat_1p_itis")
  cat(paste("iNaturalist taxa matching with ITIS",no_inat_ITIS_match," \n"))

  cat(paste("iNaturalist taxa NOT matching with ITIS",dim(inat_1p_noitis)[1]," \n"))
  parkstat$inat_ITIS_m  <- as.numeric(no_inat_ITIS_match)
  parkstat$inat_ITIS_nom  <- dim(inat_1p_noitis)[1]
  parkstat$inat_ITIS_nom_dat  <- inat_1p_noitis

  # NPS unique scientific names
  nps_1p_taxa <- as.data.frame(unique(nps_1p$ScientificName))
  nps_1p_taxa <- nps_1p[,c("ScientificName","Category")]
  names(nps_1p_taxa)[1] <- "ScientificName"
  nps_1p_taxa[] <- lapply(nps_1p_taxa, as.character)


  # Match NPS names with ITIS = 3
  nps_1p_itis <- as.data.frame(gnr_resolve(names=nps_1p_taxa$ScientificName, data_source_ids=3))
  nps_1p_itis$Category <- nps_1p$Category[match(nps_1p_itis$user_supplied_name, nps_1p$ScientificName)]
  no_nps_ITIS_match <- sqldf("select count(distinct user_supplied_name) from nps_1p_itis")
  cat(paste("NPS taxa matching with ITIS",no_nps_ITIS_match," \n"))
  nps_1p_orig$ITIS <- FALSE
  nps_1p_orig$ITIS[which(nps_1p_orig$ScientificName %in% nps_1p_itis$submitted_name)] <- TRUE

  parkstat$NPS_ITIS_m  <- as.numeric(no_nps_ITIS_match)
  nps_1p_noitis <- sqldf("select * from nps_1p_taxa where ScientificName not in (select user_supplied_name from nps_1p_itis)")
  #nps_1p_noitis$Category <- nps_1p$Category[match(nps_1p_noitis$user_supplied_name, nps_1p$ScientificName)]
  cat(paste("NPS taxa NOT matching with ITIS",dim(nps_1p_noitis)[1]," \n"))
  parkstat$NPS_ITIS_nom  <- dim(nps_1p_noitis)[1]
  parkstat$NPS_ITIS_nom_dat  <- nps_1p_noitis

  # iNaturalist taxa not in NPS

  inat_no_nps_1p <- sqldf("select * from inat_1p_taxa where taxon_species_name not in (select ScientificName from nps_1p) ")
  names(inat_no_nps_1p)[1] <- "Scientific_name"
  # Module to check synomyms in NPS
  inat_no_nps_1p_tmp <- inat_no_nps_1p[which(!inat_no_nps_1p$Scientific_name %in% syn_checked$Scientific_name),]
  if(dim(inat_no_nps_1p_tmp)[1]>0){
    inat_no_nps_1p_syn <- get_itis_binom_syn(inat_no_nps_1p$Scientific_name)
    tmpnames <- as.data.frame(inat_no_nps_1p$Scientific_name)
    names(tmpnames) <-  c("Scientific_name")
    syn_checked <- rbind(syn_checked,tmpnames)
    write.csv(syn_checked,"syn_checked.csv",row.names = F)
    inatsyn <- rbind(inatsyn,inat_no_nps_1p_syn)
    write.csv(inatsyn,"inat_syn.csv",row.names = F)
  }

  tmp_match_syn <- as.data.frame(inat_no_nps_1p_syn[which(inatsyn$Syn %in% nps_1p$ScientificName),c("Acc_name")])
  names(tmp_match_syn) <- c("ScientificName")
  nps_1p <- rbind.fill(nps_1p,tmp_match_syn)
  inat_no_nps_1p <- sqldf("select * from inat_1p_taxa where taxon_species_name not in (select ScientificName from nps_1p) ")
  names(inat_no_nps_1p)[1] <- "Scientific_name"


  cat(paste("iNaturalist taxa NOT in NPS",dim(inat_no_nps_1p)[1]," \n\n"))
  if(dim(inat_no_nps_1p)[1]>0){
    inat_no_nps_1p$iconic_taxon_name <- inat_1p$iconic_taxon_name[match(inat_no_nps_1p$Scientific_name, inat_1p$taxon_species_name)]
    # Get class as inat_no_nps_1p
    parkstat$inat_NPS_nom_dat <- getClassasIconic(inat_no_nps_1p)
    print(getClassasIconic(inat_no_nps_1p))
    inat_no_nps_1p_tab <- getIconicTab(parkstat$inat_NPS_nom_dat)
    print(inat_no_nps_1p_tab )
  } else {
    inat_no_nps_1p_tab <- NULL
  }


  # iNaturalist taxa in NPS
  inat_nps_1p <- sqldf("select * from inat_1p_taxa where taxon_species_name in (select ScientificName from nps_1p) ")
  names(inat_nps_1p)[1] <- "Scientific_name"
  inat_1p$nps <- FALSE
  inat_1p$nps[which(inat_1p$taxon_species_name %in% inat_nps_1p$Scientific_name)] <- TRUE
  inat_1p$NPS_Occurrence <- nps_1p[match(inat_1p$taxon_species_name, nps_1p$ScientificName),c("Occurrence")]
  inat_1p$NPS_RecordStatus <- nps_1p[match(inat_1p$taxon_species_name, nps_1p$ScientificName),c("RecordStatus")]
  inat_1p$NPS_Nativeness <- nps_1p[match(inat_1p$taxon_species_name, nps_1p$ScientificName),c("Nativeness")]
  inat_1p$NPS_Category <- nps_1p[match(inat_1p$taxon_species_name, nps_1p$ScientificName),c("Category")]
  inat_1p$NPS_CommonNames <- nps_1p[match(inat_1p$taxon_species_name, nps_1p$ScientificName),c("CommonNames")]
  inat_1p$NPS_GRank <- nps_1p[match(inat_1p$taxon_species_name, nps_1p$ScientificName),c("GRank")]
  inat_1p$focal <- FALSE
  inat_1p$focal[which(inat_1p$iconic_taxon_name %in% c("Amphibia", "Aves", "Mammalia", "Reptilia", "Plantae"))] <- TRUE

  nps_1p_orig$inat <- FALSE
  nps_1p_orig$inat[which(nps_1p_orig$ScientificName %in% inat_nps_1p$Scientific_name)] <- TRUE

  cat(paste("iNaturalist taxa in NPS",dim(inat_nps_1p)[1]," \n\n"))
  inat_nps_1p$iconic_taxon_name <- nps_1p$Category[match(inat_nps_1p$Scientific_name,
                                                         nps_1p$ScientificName)]
  iconic_taxon_name <- NULL
  Scientific_name <- NULL
  inat_nps_1p_tab <- as.data.frame(dplyr::summarise(dplyr::group_by(inat_nps_1p,
                                                                    iconic_taxon_name),
                                                    ct=dplyr::n_distinct(Scientific_name)))
  print(inat_nps_1p_tab )

  parkstat$inat_NPS_nom  <- dim(inat_no_nps_1p)[1]
  parkstat$inat_NPS_nom_sum   <- inat_no_nps_1p_tab
  parkstat$inat_NPS_nom_dat <- inat_no_nps_1p
  parkstat$inat_NPS_match  <- dim(inat_nps_1p)[1]
  parkstat$inat_NPS_match_sum <- inat_nps_1p_tab
  parkstat$inat <- inat_1p
  parkstat$NPS <- nps_1p_orig
  return(parkstat)
}
