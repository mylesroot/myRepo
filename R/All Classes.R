#' The PlotProfile object
#'
#' @export
#'
#'

setClass("ChIPprofile",contains = "RangedSummarizedExperiment",
         slots=c(params="list"
         ))


#' Edited Cbind
#' @rdname manipulateObjects
#' @export
setMethod("cbind", "ChIPprofile",
          function (x,...,deparse.level=1)
          {
            assayList <- vector("list",length=length(assays(x)))
            regions <- list(x,...)
            for(a in 1:length(assays(x))){
              listTemp <- vector("list",length=length(regions))
              for(r in 1:length(regions)){
                listTemp[[r]] <- assays(regions[[r]])[[a]]
              }
              assayList[[a]] <- do.call(cbind,listTemp)
            }
            subsetProfile <- SummarizedExperiment(assayList,rowRanges=rowRanges(x))
            metadata(subsetProfile)$names <- metadata(x)$names
            metadata(subsetProfile)$AlignedReadsInBam <- metadata(x)$AlignedReadsInBam
            return(new("ChIPprofile",subsetProfile,params=x@params))
          }
)
