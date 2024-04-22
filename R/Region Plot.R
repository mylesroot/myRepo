library (rtracklayer)
library(genomation)
library(soGGi)

bw1forward = ("~/Documents/Shaun Data/CRAC7_4th_repeat_NNNGCGCAGC_L5Ac_collapsed.bs1.forward.bw")
bw1reverse = ("~/Documents/Shaun Data/CRAC7_4th_repeat_NNNGCGCAGC_L5Ac_collapsed.bs1.reverse.bw")
bw2forward = ("~/Documents/Shaun Data/CRAC7_4th_repeat_NNNATTAGC_L5Ab_collapsed.bs1.forward.bw")
bw2reverse = ("~/Documents/Shaun Data/CRAC7_4th_repeat_NNNATTAGC_L5Ab_collapsed.bs1.reverse.bw")
bw3forward = ("~/Documents/Shaun Data/CRAC7_4th_repeat_NNNTAAGC_L5Aa_collapsed.bs1.forward.bw")
bw3reverse = ("~/Documents/Shaun Data/CRAC7_4th_repeat_NNNTAAGC_L5Aa_collapsed.bs1.reverse.bw")

# L5ab Sample Name = "CRAC7_4_RNH201-HTP"

# L5ac Sample Name = CRAC7_4_RPO21-HTP

# L5aa Sample Name = "CRAC7_4_untagged"


# The function and it's variables

# I'm assuming inputs are always going to be a bigWig file and a GRanges object and that the method is always going to be equal bins.

regionPlotEdited = function(bigWig, testRanges, nOfWindows, samplename=NULL, feature = "tss", upstreamBinNumber = 10, regionBinNumber = 50, style = "percentOfRegion") {





# Input is bigWig and GRanges (remove all other input possibilites)

# SoGGi then checks input parameters. If the style is percentOfregion then it sets distanceAround to 100 and distanceUp and distanceDown to distanceAround. It sets distanceAround to 1500 as well as distanceUp and distanceDown to the value of
# distance around if these values are NULL.





## Initialize empty matrices and paramaters for collecting coverage analysis
## Find maximum distance to use for filtering out of bounds extended GRanges

# if(style == "percentOfRegion"){

  # Edited code here as I don't need the distance around variable and will just be calculating for regions within the file
  maxDistance = width(testRanges)
  # maxDistance <- round((distanceAround/100)*width(testRanges))

  RegionsMat <- NULL
  distanceUpStart <- NULL
  distanceDownEnd <- NULL

totalReads <- NA




# Import data and find contig information


#if(format=="bigwig"){
  message("Reading BigWig contig information...",appendLF = FALSE)
  bwFF <- BigWigFile(bigWig)
  lengths <- seqlengths(bwFF)
  allchrs <- names(lengths)
  message("..Done")
#}



# Exclude and count regions beyond contig boundaries

  # if(style == "percentOfRegion"){
    message("Filtering regions which extend outside of genome boundaries...",appendLF = FALSE)
    testRangeNames <- unique(seqnames(testRanges))
    temptestranges <- GRanges()
    for(i in 1:length(testRangeNames)){
      perChrMaxDistance <- maxDistance[as.vector(seqnames(testRanges) %in% as.vector(testRangeNames[i]))]
      perchrRanges <- testRanges[seqnames(testRanges) %in% as.vector(testRangeNames[i])]
      temptestranges <- c(temptestranges,perchrRanges[end(perchrRanges)+perChrMaxDistance < lengths[names(lengths) %in% testRangeNames[i]]
                                                      & start(perchrRanges)-perChrMaxDistance > 0 ])
      # print(i)
      perChrMaxDistance <- perChrMaxDistance[end(perchrRanges)+perChrMaxDistance < lengths[names(lengths) %in% testRangeNames[i]]
                                             & start(perchrRanges)-perChrMaxDistance > 0 ]
      distanceUpStart <- c(distanceUpStart,perChrMaxDistance)
    }
    distanceDownEnd <- distanceUpStart

  #}
    message("..Done")
    message("Filtered ",length(testRanges)-length(temptestranges)," of ",length(testRanges)," regions")
    testRanges <- temptestranges
    temptestranges <- NULL






# Split ranges into +/- strand. Regions with no strand information are assigned to + strand




message("Splitting regions by Watson and Crick strand..",appendLF = FALSE)

# Add a specific giID column to each range of the testRanges object

mcols(testRanges) <- cbind(mcols(testRanges),data.frame(giID = paste0("giID",seq(1,length(testRanges)))))

# Assign unknown strands to positive and split the testRanges object into two different GRanges objects with + or - strands.

strand(testRanges[strand(testRanges) == "*"]) <- "+"
testRangesPos <- testRanges[strand(testRanges) == "+"]
testRangesNeg <- testRanges[strand(testRanges) == "-"]
message("..Done")

# Code below removed as there's no style variable for my function

# if(style=="percentOfRegion"){
#  distanceUpStartPos <- distanceUpStart[as.vector(strand(testRanges) == "+")]
#  distanceDownEndPos <- distanceUpStartPos
#  distanceUpStartNeg <- distanceUpStart[as.vector(strand(testRanges) == "-")]
#  distanceDownEndNeg <- distanceUpStartNeg
#  message("..Done")
#}else{
#  distanceUpStartPos <- distanceUpStart
#  distanceDownEndPos <- distanceDownEnd
#  distanceUpStartNeg <- distanceUpStart
#  distanceDownEndNeg <- distanceDownEnd
#  message("..Done")
#}

message("Found ",length(testRangesPos)," Watson strand regions")
message("Found ",length(testRangesNeg)," Crick strand regions")

## Extend regions and get positive versus negative regions
# I will only be getting pos vs negative

message("Concatenating pos and neg Ranges..",appendLF=FALSE)

exttestRanges = c(testRangesPos, testRangesNeg)

# exttestRanges <- c(GRanges(seqnames(testRangesPos),IRanges(start(testRangesPos)-distanceUpStartPos,end(testRangesPos)+distanceDownEndPos)),
#                  GRanges(seqnames(testRangesNeg),IRanges(start(testRangesNeg)-distanceDownEndNeg,end(testRangesNeg)+distanceUpStartNeg))
#)
message("...done")


## Create GRanges to be used in scanBamParam while reading in Bamfile regions.
reducedExtTestRanges <- reduce(exttestRanges)

# if(format=="bigwig"){
  message("Extracting coverage info from BigWig file and importing as RleList..",appendLF = FALSE)
  bwSelect <- BigWigSelection(GRanges(seqnames=seqnames(reducedExtTestRanges[seqnames(reducedExtTestRanges) %in% allchrs]),IRanges(start=start(reducedExtTestRanges[seqnames(reducedExtTestRanges) %in% allchrs]),end=end(reducedExtTestRanges[seqnames(reducedExtTestRanges) %in% allchrs]))))
  genomeCov <- import.bw(bigWig,selection=bwSelect,as="RleList")
# }


# Filters chromosomes to include only those chromosomes present within the GRanges object named reducedExtTestRanges.

chromosomes <- seqlevels(genomeCov)
chromosomes <- chromosomes[chromosomes %in% unique(seqnames(reducedExtTestRanges))]


# Next steps assume that the method is Percent of Region and that we're using the Bin as opposed to Spline argument

meansListNeg <- vector("numeric")
meansListPos <- vector("numeric")

grListWindowsPos <- GRanges()
grListWindowsNeg <- GRanges()

## Create GRanges of windows across regions
message("Making windows.")

## Positive regions

if(length(testRangesPos) > 0){

## Calculate bin lengths
grWidths <- width(testRangesPos)
windows <- floor(grWidths%/%nOfWindows)
extraForWindows <- grWidths%%nOfWindows
# extraForFlankWindows <- grWidths%%(nOfWindows*((distanceAround)/100))
addToWindow <- 0
startPos <- start(testRangesPos)#-distanceUpStartPos#-(windows/2)

# Creates a vector rem filled with zeros, the same length as extraForFlankWindows. Not needed I don't think.

# rem <- rep(0,length(extraForFlankWindows))
# rem2 <- NULL



## Create bin GRanges for positive 5' flanking regions. NOT NEEDED - HERE FOR REFERENCE.

#message("Windowing positive 5' flanking ")
#for(i in 1:(nOfWindows*((distanceAround)/100))){
  #message("Window", i,appendLF=F)
#  rem2 <- rem+((extraForFlankWindows >= i)+0)

#  grListWindowsPos <- c(grListWindowsPos,GRanges(seqnames(testRangesPos),IRanges(
#    (startPos)+rem+(windows*(i-1)),
#    startPos+(windows*i)-1+rem2),giID=testRangesPos$giID))

#  rem <- rem2


## Create bin GRanges for positive regions

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

## Create bin GRanges for positive 3' flanking regions - NOT NEEDED, HERE FOR REFERENCE

#rem <- rep(0,length(extraForFlankWindows))
#rem2 <- NULL
#startPos <- end(testRangesPos)#-(windows/2)
#message("Windowing positive 3' flank ")
#for(i in 1:(nOfWindows*((distanceAround)/100))){
  #message("Window", i)
#  rem2 <- rem+((extraForFlankWindows >= i)+0)

#  grListWindowsPos <- c(grListWindowsPos,GRanges(seqnames(testRangesPos),IRanges(
#    (startPos)+rem+(windows*(i-1)),
#    startPos+(windows*i)-1+rem2),giID=testRangesPos$giID))
#  rem <- rem2
#}

## Order by giID to group windows from gene. Retains secondary order as
## created (bin order)
grListWindowsPos <- grListWindowsPos[order(grListWindowsPos$giID)]

}


## Handle negative GRanges as with positive GRanges
if(length(testRangesNeg) > 0){

  grWidths <- width(testRangesNeg)
  windows <- floor(grWidths%/%nOfWindows)
  extraForWindows <- grWidths%%nOfWindows
#  extraForFlankWindows <- grWidths%%(nOfWindows*((distanceAround)/100))
  addToWindow <- 0
  startPos <- start(testRangesNeg)#-distanceDownEndNeg
#  rem <- rep(0,length(extraForFlankWindows))
  rem2 <- NULL


#  message("Windowing negative 5' flank ")
#  for(i in 1:(nOfWindows*((distanceAround)/100))){
    #message("Window", i)
#    rem2 <- rem+((extraForFlankWindows >= i)+0)

#    grListWindowsNeg <- c(grListWindowsNeg,GRanges(seqnames(testRangesNeg),IRanges(
#      (startPos)+rem+(windows*(i-1)),
#      startPos+(windows*i)-1+rem2),giID=testRangesNeg$giID))
#    rem <- rem2
#}

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


#  message("Windowing negative 3' flank ")

#  for(i in 1:(nOfWindows*((distanceAround)/100))){
    #message("Window", i)
#    rem2 <- rem+((extraForFlankWindows >= i)+0)

#    grListWindowsNeg <- c(grListWindowsNeg,GRanges(seqnames(testRangesNeg),IRanges(
#      (startPos)+rem+(windows*(i-1)),
#      startPos+(windows*i)-1+rem2),giID=testRangesNeg$giID))
#    rem <- rem2
#  }

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


# Stuff after this stage is likely where problems arise (I don't know what the code below does exactly)

## Create matrices for mean bin coverage

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
#print(profileMat)



### Create ChIPprofile object.

colnames(profileMat) <- c(paste0(feature, " Bin ", seq(1,nOfWindows)))
filteredRanges <- c(testRangesPos,testRangesNeg)


# This alters the data that's returned to the assay part of the ChIProfile
# Commented out code is SoGGi and works

# profileSample <- SummarizedExperiment(profileMat,rowRanges=filteredRanges[match(rownames(profileMat),filteredRanges$giID)])



profileSample = SummarizedExperiment(profileMat, rowRanges=filteredRanges[match(rownames(profileMat), filteredRanges$giID)])




# Not sure what this code does - leaving here as reference

# if(is.null(samplename)){
#  if(format %in% c("rlelist","pwm","granges")){
#    metadata(profileSample)  <- list(names=c("Sample"))
#  }else{
#    metadata(profileSample)<- list(names=c(bamFile),AlignedReadsInBam=totalReads)
#  }
#}else{

  metadata(profileSample)<- list(names=samplename,AlignedReadsInBam=totalReads)
#}

## Pass parameters

paramList <- list("nOfWindows"=nOfWindows,
                  "style"=style, "upstreamBinNumber" = upstreamBinNumber, "regionBinNumber" = regionBinNumber)
#                  "samplename"=samplename,
#                  "nOfWindows"=nOfWindows,
#                  "FragmentLength"=FragmentLength,
#                  "distanceAround"=distanceAround,
#                  "distanceUp"=distanceUp,
#                  "distanceDown"=distanceDown,
#                  "distanceInRegionStart"=distanceInRegionStart,
#                  "distanceInRegionEnd"=distanceInRegionEnd,
#                  "distanceOutRegionStart"=distanceOutRegionStart,
#                  "distanceOutRegionEnd"=distanceOutRegionEnd,
#                  "paired"=paired,
#                  "normalize"=normalize,
#                  "plotBy"=plotBy,
#                  "removeDup"=removeDup,
#                  "format"=format,
#                  "seqlengths"=seqlengths,
#                  "forceFragment"=forceFragment,
#                  "method"=method,
#                  "genome"=genome,
#                  "cutoff"=cutoff,
#                  "minFragmentLength"=minFragmentLength,
#                  "maxFragmentLength"=maxFragmentLength,
#                  "downSample"=downSample
#)
return(new("ChIPprofile",profileSample,params=paramList))

}
