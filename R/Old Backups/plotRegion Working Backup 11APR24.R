install.packages("ggplot2")
library(ggplot2)

install.packages("reshape2")
library(reshape2)


# @param gts A list of character vectors or GRangesList

plotRegionEdited = function(object, summariseBy = c("name")){

  # Pass parameters for bin sizes (used for break calculation)

  upstreamBinNumber = object@params$upstreamBinNumber
  regionBinNumber = object@params$regionBinNumber



  # Creates profileList - creates a list of the different assays (in matrix-like format) in the ChiProfile object. For each assay, it calculates the mean
  #                       score for the columns (bins)


  profileList <- lapply(c(assays(object)),colMeans,na.rm=T)

  print(profileList)
  print("1")
  #  }

  # Creates profileFrame. If I have two or more assays in profileList, this combines them side-by-side into a single matrix.
  # For the column names of each matrix, it unlists the "names" metadata list of the ChIProfile object. The "names" metadata list is
  # created from the samplename parameter in the regionPlot function.

  profileFrame <- do.call(cbind,profileList)
  colnames(profileFrame) <- basename(unlist(metadata(object)["names"]))
  print(profileFrame)
  print("2")

  ## Attach index for different styles of plots

  # SOGGI IS USING PARAMETER INFORMATION AT THIS STAGE.

  # Creates axisIndex which is a list of integers from 1-length of the profile frame (i.e. how many bins there are).

  # axisIndex=c(seq(1,object@params$nOfWindows))
  axisIndex=c(seq(1,nrow(profileFrame)))
  # nOfWindows = object@params$nOfWindows
  nOfWindows = nrow(profileFrame)
  print(axisIndex)
  print("3")

  ## Reformat profileFrame data into long format. The rowNames become the axisIndex values. The colNames become xIndex, Sample, and score.

  rownames(profileFrame) <- axisIndex
  print(rownames(profileFrame))
  print("4")
  meltedProfileFrame <- melt(profileFrame)
  colnames(meltedProfileFrame) <- c("xIndex","Sample","Score")
  print(meltedProfileFrame)


  P <- ggplot(meltedProfileFrame,
              aes_string(x="xIndex",y="Score"))+geom_path(alpha = 1,size=1.3)+xlim(0,max(axisIndex))+ylab("Score")+theme(axis.title.y=element_text(angle=0))

  labels = rownames(profileFrame)


  P <- P + scale_x_continuous(breaks=c(1, upstreamBinNumber, (regionBinNumber+upstreamBinNumber), (regionBinNumber+regionBinNumber+upstreamBinNumber), (regionBinNumber+regionBinNumber+regionBinNumber+upstreamBinNumber), (regionBinNumber+regionBinNumber+regionBinNumber+regionBinNumber+upstreamBinNumber)),
                              labels=c("Start-", "TSS", "End Of Exon 1", "End Of Intron 1", "TES", "End-")
  )+
    theme(axis.text.x  = element_text(angle=0, vjust=0.5, size=6))


  #  P <- P + scale_x_continuous(breaks=c(nOfWindows),
  #                              labels=c(paste0("Tss")
  #                              ))+
  #    theme(axis.text.x  = element_text(angle=45, vjust=0.5, size=12))


  P <- P+aes_string(colour="Sample")


  return(P)

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

