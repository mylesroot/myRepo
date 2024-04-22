#' Function To Plot Your Profile Plot objects
#'
#' This function takes the Profile Plot objects you've created and enables you to plot them.
#'
#' @param object Profile Plot object you created
#' @param group = A list of gene names you'd like to facet your plots with
#' @param goupBy = A character vector of the group names to display on your plot
#'
#' @import BiocGenerics reshape2 genomation GenomicRanges ggplot2
#'
#' @export


plotProfile = function(object, group = NULL, groupBy = NULL){

  # Pass parameters for bin sizes (used for break calculation)

  breaks = cumsum(object@params$windowLengths)
  maxBreaks = length(object@params$windowLengths)
  levels = list()


  # Creates profileList - creates a list of the different assays (in matrix-like format) in the ChiProfile object. For each assay, it calculates the mean
  #                       score for the columns (bins)

  if(!is.null(group)) {

    profileList <- list()

    for(p in 1:length(assays(object))){

      ## Create a list of profile matrices
      profileTemp <- assays(object)[[p]]

      profileListTemp = list()

      found_genes <- rowRanges(object)$name[rowRanges(object)$name %in% group]

      found_gene_IDs = rowRanges(object)$giID[rowRanges(object)$name %in% found_genes]

      split_index = rownames(profileTemp) %in% found_gene_IDs

      group1 = profileTemp[split_index, ]

      group2 = profileTemp[!split_index,]

      profileListTemp = c(profileListTemp, list(group1, group2))

      profileListTemp <- lapply(profileListTemp,colMeans,na.rm=T)

      profileMatTemp <- melt(as.data.frame(do.call(cbind,profileListTemp)))

      profileMatTemp$variable = factor(profileMatTemp$variable, levels = c("V1", "V2"), labels = groupBy)

      #axisIndex=c(seq(1,nrow(profileMatTemp)))

      axisIndex = c(seq(1,sum(object@params$windowLengths)))

      profileFrame <-data.frame("xIndex"=axisIndex,Group=profileMatTemp[,1],Sample=basename(unlist(metadata(object)["names"]))[p],Score=profileMatTemp[,2])

      profileList[[p]] <- profileFrame

    }

    meltedProfileFrame <- do.call(rbind,profileList)

    colnames(meltedProfileFrame) <- c("xIndex","Group","Sample","Score")


  } else {

    profileList <- lapply(c(assays(object)),colMeans,na.rm=T)

    profileFrame <- do.call(cbind,profileList)

    colnames(profileFrame) <- basename(unlist(metadata(object)["names"]))

    axisIndex=c(seq(1,nrow(profileFrame)))

    nOfWindows = nrow(profileFrame)

    rownames(profileFrame) <- axisIndex

    meltedProfileFrame <- melt(profileFrame)

    colnames(meltedProfileFrame) <- c("xIndex","Sample","Score")

  }

  P <- ggplot(meltedProfileFrame,
              aes_string(x="xIndex",y="Score"))+geom_path(alpha = 1,size=1.3)+xlim(0,max(axisIndex))+ylab("Score")+theme(axis.title.y=element_text(angle=0))

  labels = rownames(profileFrame)


  P <- P + scale_x_continuous(breaks=c(1, as.numeric(paste0(breaks))),
                              labels=c("-2000", "TSS", "5'SS", "3'SS", "PolyA", "+2000")#lapply(breaks, function(x) paste0("Break")))
  )+
    theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=10))


  #  P <- P + scale_x_continuous(breaks=c(nOfWindows),
  #                              labels=c(paste0("Tss")
  #                              ))+
  #    theme(axis.text.x  = element_text(angle=45, vjust=0.5, size=12))


  P <- P+aes_string(colour="Sample")

  if (!is.null(group)){

    #    if (!is.null(groupBy)){

    #    facet <- facet_wrap(
    #      formula(paste("~",paste(groupBy,collapse="+"))), scales = "free_y"
    #    )

    #    }
    #    else {

    facet <- facet_wrap(
      formula(paste("~",paste("Group",collapse="+"))), scales = "free_y"
    )


    P <- P + facet

  }

  return(P)

}


