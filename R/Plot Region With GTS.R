rpGenes = "~/Documents/Shaun Data/RP_Intron_Genes.bed"
rpGenes = readBed(rpGenes)

# Get the names in rpGenes
# Find out what the giIDs of these names are in my assay
# Check if the giIDs match any in rownames(assay(object))
# If they do, create a separate profile for them

names = rpGenes$name
gts = list()

found_genes <- rowRanges(object)$name[rowRanges(object)$name %in% names]

found_gene_IDs = rowRanges(test)$giID[rowRanges(test)$name %in% found_genes]

gts = found_gene_IDs


plotRegionTesting <- function(object, gts=gts, summariseBy=NULL)
{


  ## When running with gts option
  ## SummariseBy can now only be used to select column for gts to be match to
  ## and colourBy, lineBy and groupBy will refer to sample metadata or be "group"

  if(!is.null(gts)){

    profileList <- list()

    ## Start cycling through assays
    for(p in 1:length(assays(object))){

      ## extract profile matrix
      profileTemp <- assays(object)[[p]]

      profileTempList <- lapply(gts,function(x)colMeans(subsetProfile(profileTemp,x),na.rm=T))

    }
  }
}



subsetProfile = function(profile, group){

  return(profile[rownames(profile) %in% group,])

}

subsetProfile <- function(profile,group,granges,summariseColumn){
    if(class(group) == "character"){
      if(is.null(summariseColumn)){
        return(profile[rownames(profile) %in% group,])
      }else{
        return(profile[mcols(granges)[,summariseColumn] %in% group,])
      }
    }
    if(class(group) == "GRanges"){
      return(profile[granges %over% group,])
    }
  }









