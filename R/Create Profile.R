#' Create a summarised experiment object from a list of GRanges and a .bigWig file.
#'
#' This function takes a bigWig, GRangesList, list of window lengths, and a sample name.
#' It returns a summarised experiment object that can then be used for plotting.
#'
#'@param bigWig A file path to a .bigWig file
#'@param testRangesList A GRangesList object containing the gene ranges to calculate scores over
#'@param nOfWindowsList A numeric vector containing the window lengths for binning. Each length in the list
#'must correspond to the same position in the GRanges list.
#'@param samplename The name of your sample, used to group data in the plot.
#'
#'@import reshape2 BiocGenerics GenomicRanges genomation rtracklayer SummarizedExperiment
#'
#'@export


createProfile = function(bigWig, testRangesList, nOfWindowsList, samplename=NULL) {

  chipProfiles = NULL

  # Check lists are the same length

  if (length(testRangesList) != length(nOfWindowsList)) {
    stop("List of window lengths should match the length of the testRangesList.")
  }

  for(i in seq_along(testRangesList)){

    testRanges = testRangesList[[i]]

    nOfWindows = nOfWindowsList[[i]]

    maxDistance = width(testRanges)

    RegionsMat <- NULL
    distanceUpStart <- NULL
    distanceDownEnd <- NULL

    totalReads <- NA

    message("Reading BigWig contig information...",appendLF = FALSE)
    bwFF <- BigWigFile(bigWig)
    lengths <- seqlengths(bwFF)
    allchrs <- names(lengths)
    message("..Done")

    message("Splitting regions by Watson and Crick strand..",appendLF = FALSE)

    mcols(testRanges) <- cbind(mcols(testRanges),data.frame(giID = paste0("giID",seq(1,length(testRanges)))))

    strand(testRanges[strand(testRanges) == "*"]) <- "+"
    testRangesPos <- testRanges[strand(testRanges) == "+"]
    testRangesNeg <- testRanges[strand(testRanges) == "-"]
    message("..Done")

    message("Found ",length(testRangesPos)," Watson strand regions")
    message("Found ",length(testRangesNeg)," Crick strand regions")

    message("Concatenating pos and neg Ranges..",appendLF=FALSE)

    exttestRanges = c(testRangesPos, testRangesNeg)

    message("...done")

    reducedExtTestRanges <- reduce(exttestRanges)

    message("Extracting coverage info from BigWig file and importing as RleList..",appendLF = FALSE)
    bwSelect <- BigWigSelection(GRanges(seqnames=seqnames(reducedExtTestRanges[seqnames(reducedExtTestRanges) %in% allchrs]),IRanges(start=start(reducedExtTestRanges[seqnames(reducedExtTestRanges) %in% allchrs]),end=end(reducedExtTestRanges[seqnames(reducedExtTestRanges) %in% allchrs]))))
    genomeCov <- import.bw(bigWig,selection=bwSelect,as="RleList")

    chromosomes <- seqlevels(genomeCov)
    chromosomes <- chromosomes[chromosomes %in% unique(seqnames(reducedExtTestRanges))]

    meansListNeg <- vector("numeric")
    meansListPos <- vector("numeric")

    grListWindowsPos <- GRanges()
    grListWindowsNeg <- GRanges()

    message("Making windows.")

    if(length(testRangesPos) > 0){

      grWidths <- width(testRangesPos)
      windows <- floor(grWidths%/%nOfWindows)
      extraForWindows <- grWidths%%nOfWindows
      addToWindow <- 0
      startPos <- start(testRangesPos)

      startPos <- start(testRangesPos)#-(windows/2)
      rem <- rep(0,length(extraForWindows))
      rem2 <- NULL
      message("Windowing positive regions ")
      for(i in 1:(nOfWindows)){
        #message("Window", i)
        rem2 <- rem+((extraForWindows >= i)+0)

        grListWindowsPos <- c(grListWindowsPos,GRanges(seqnames(testRangesPos),IRanges(
          (startPos)+rem+(windows*(i-1)),
          startPos+(windows*i)-1+rem2),giID=testRangesPos$giID))
        rem <- rem2
      }

      grListWindowsPos <- grListWindowsPos[order(grListWindowsPos$giID)]

    }

    if(length(testRangesNeg) > 0){

      grWidths <- width(testRangesNeg)
      windows <- floor(grWidths%/%nOfWindows)
      extraForWindows <- grWidths%%nOfWindows
      #  extraForFlankWindows <- grWidths%%(nOfWindows*((distanceAround)/100))
      addToWindow <- 0
      startPos <- start(testRangesNeg)#-distanceDownEndNeg
      #  rem <- rep(0,length(extraForFlankWindows))
      rem2 <- NULL


      startPos <- start(testRangesNeg)#-(windows/2)
      rem <- rep(0,length(extraForWindows))
      rem2 <- NULL
      message("Windowing negative regions ")

      for(i in 1:(nOfWindows)){
        #message("Window", i)
        rem2 <- rem+((extraForWindows >= i)+0)

        grListWindowsNeg <- c(grListWindowsNeg,GRanges(seqnames(testRangesNeg),IRanges(
          (startPos)+rem+(windows*(i-1)),
          startPos+(windows*i)-1+rem2),giID=testRangesNeg$giID))
        rem <- rem2
      }

      #  rem <- rep(0,length(extraForFlankWindows))
      rem2 <- NULL
      startPos <- end(testRangesNeg)#-(windows/2)


      grListWindowsNeg <- grListWindowsNeg[order(grListWindowsNeg$giID)]

    } else {

      message("No negative regions found")

    }

    grListWindows <- list(grListWindowsPos,grListWindowsNeg)
    message("..done\n")


    ## Cycle through contigs to extract scores from rlelist per contig

    message(paste0("Calculating bin scores for regions.\nProcessing per contig"))


    for(c in 1:length(chromosomes)){

      message(paste0("contig: ",c))

      message("Processing inner region windows in ",chromosomes[c])
      covPerPeakPos <- Views(genomeCov[[which(names(genomeCov) %in% chromosomes[c])]],ranges(grListWindows[[1]][seqnames(grListWindows[[1]]) == chromosomes[c]]))
      doubleTempPos <- viewMeans(covPerPeakPos)
      names(doubleTempPos) <- as.vector(grListWindows[[1]][seqnames(grListWindows[[1]]) == chromosomes[c]]$giID)
      meansListPos <- c(meansListPos,doubleTempPos)
      covPerPeakNeg <- Views(genomeCov[[which(names(genomeCov) %in% chromosomes[c])]],ranges(grListWindows[[2]][seqnames(grListWindows[[2]]) == chromosomes[c]]))
      doubleTempNeg <- viewMeans(covPerPeakNeg)
      names(doubleTempNeg) <- as.vector(grListWindows[[2]][seqnames(grListWindows[[2]]) == chromosomes[c]]$giID)
      meansListNeg <- c(meansListNeg,doubleTempNeg)
      message("..done")

    }


    meansPos <- matrix(meansListPos,
                       ncol=nOfWindows
                       ,byrow=TRUE)
    if(nrow(meansPos) > 0){
      rownames(meansPos) <- matrix(names(meansListPos),ncol=nOfWindows
                                   ,byrow=TRUE)[,1]
    }
    meansNeg <- matrix(meansListNeg,
                       ncol=nOfWindows
                       ,byrow=TRUE)[,nOfWindows:1]
    if(nrow(meansNeg) > 0){

      rownames(meansNeg) <- matrix(names(meansListNeg),ncol=nOfWindows
                                   ,byrow=TRUE)[,1]
    }
    meansMat <- rbind(meansPos,meansNeg)
    profileMat <- meansMat[order(rownames(meansMat)),]


    colnames(profileMat) <- c(paste0("Bin ", seq(1,nOfWindows)))
    filteredRanges <- c(testRangesPos,testRangesNeg)

    profileSample = SummarizedExperiment(profileMat, rowRanges=filteredRanges[match(rownames(profileMat), filteredRanges$giID)])

    metadata(profileSample)<- list(names=samplename,AlignedReadsInBam=totalReads)

    paramList <- list("nOfWindows"=nOfWindows, "windowLengths" = nOfWindowsList)


    newChip = (new("ChIPprofile",profileSample,params=paramList))
    chipProfiles = c(chipProfiles, newChip)

  }

  joinedProfiles = do.call(cbind, chipProfiles)
  return(joinedProfiles)
  #return(new("ChIPprofile",profileSample,params=paramList))
}





