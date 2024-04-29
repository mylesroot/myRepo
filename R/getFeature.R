#' Get a specific feature of interest
#'
#' Input the strand, start and end point of the range, and specify the feature you want this around. Speficy upstream
#' and downstream as well.
#'
#'@param object The GRanges object to define your ranges from.
#'@param start_feature The feature you would like to start around.
#'@param start_flank How many bp to set as your flank region around the start.
#'@param start_exon If your start feature is an exon/intron, set this number to whichever you'd like to start around.
#'@param start_exon_boundary The exon boundary you'd like to start from. Set to '3prime' or '5prime'
#'@param start_direction The direction from the start point to set your flank region start.
#'@param end_feature The feature you'd like to end at.
#'@param end_flank How many bp to set as your flank region around the end point.
#'@param end_direction The direction from the end point to set your flank region start.
#'@param end_exon If your end feature is an exon/intron, set this number to whichever you'd like to end around.
#'@param end_exon_boundary The exon boundary you'd like to end at. Set to '3prime' or '5prime'
#'@param return Set the strandedness of the ranges you'd like to return. This can be set to 'all', 'plus', or 'minus'.
#'
#'@importFrom GenomicRanges GRanges
#'@importFrom GenomicRanges GRangesList
#'@import BiocGenerics
#'
#'@export

#Possible features: "TSS", "TES", "Exon"

getFeature = function(object, start_feature = "TSS", start_flank = 0, start_exon = NULL, start_exon_boundary = NULL,
                      start_direction = "up", end_feature = "TES", end_flank = 0, end_direction = 'down',
                      end_exon = NULL, end_exon_boundary = NULL, return = 'all'){

  # Ensure input validity
  if (!inherits(object, "GRanges")) {
    stop("Input must be a GRanges object.")
  }

  allowedFeatures = c("TSS", "TES", "Exon")

  if(!start_feature %in% allowedFeatures){

    stop("start_feature must be either 'TSS', 'TES' or 'Exon'.")

  }

  if(!end_feature %in% allowedFeatures){

    stop("end_feature must be either 'TSS', 'TES' or 'Exon'.")

  }

  allowed_directions = c("up", "down")

  if(!start_direction %in% allowed_directions){

    stop("start_direction should be set to either 'up' or 'down'.")

  }

  if(!end_direction %in% allowed_directions){

    stop("stop_direction should be set to either 'up' or 'down'")

  }

  # Define empty variables

  starts = list()
  object = object
  plusStrandRanges = GRanges()
  minusStrandRanges = GRanges()

  # If Exon is The Start, first filter and then separate the strands

  if(start_feature == "Exon"){

    # First check that the boundary is not set to NULL

    if(!is.null(start_exon_boundary)){

      # Get the exon number

      exonStartNumber = start_exon

      # Check to see if there are genes with that many exons

      splitIndex = object$blockCount >= exonStartNumber

      # If there are, filter the bed only to include those ones.

      if(TRUE %in% splitIndex){

        object = object[splitIndex]

        # Now check that the end feature is not an exon - if it is we filter

        if(end_feature != "Exon"){

          # split the strands in object for processing

          plusStrandRanges <- subset(object, strand(object) == "+")

          minusStrandRanges <- subset(object, strand(object) == "-")

        } else {

          if(!is.null(end_exon)){

            exonEndNumber = end_exon

            splitIndex = object$blockCount >= exonEndNumber

            if(TRUE %in% splitIndex){

              object = object[splitIndex]

              plusStrandRanges <- subset(object, strand(object) == "+")

              minusStrandRanges <- subset(object, strand(object) == "-")

            } else {

              stop("There are no genes with that many exons. Please enter another exon end.")

            }
          } else {

            stop("If end_feature is set to Exon then the end_exon must not be null.")

          }
        }

      } else {

          stop("There are no genes with that many exons. Please enter another exon start.")

        }

    } else {

      stop("If the start feature is set to Exon, start_exon_boundary must be set to either '3prime' or '5prime'.")

    }

  }

    else if(start_feature == "TSS"){

      # Check if there is no exon start information

      if(is.null(start_exon) & is.null(start_exon_boundary)){

        if(end_feature != "Exon"){

          # split the strands in object for processing

          plusStrandRanges <- subset(object, strand(object) == "+")

          minusStrandRanges <- subset(object, strand(object) == "-")

        } else {

          if(!is.null(end_exon)){

            # Get the exon number

            exonEndNumber = end_exon

            splitIndex = object$blockCount >= exonEndNumber

            if(TRUE %in% splitIndex){

              object = object[splitIndex]

              plusStrandRanges <- subset(object, strand(object) == "+")

              minusStrandRanges <- subset(object, strand(object) == "-")

            } else {

              stop("There are no genes with that many exons. Please enter another exon end.")

            }
          } else {

            stop("If end_feature is set to Exon then end_exon must not be null")

          }
        }

      } else {

        stop("If start feature is TSS then no exon start information must be entered.")

      }

    }



  # If it's TES, set start to TES

  else {

    if(is.null(start_exon) & is.null(start_exon_boundary)){

      if(end_feature != "Exon"){

        # split the strands in object for processing

        plusStrandRanges <- subset(object, strand(object) == "+")

        minusStrandRanges <- subset(object, strand(object) == "-")

      } else {

        if(!is.null(end_exon)){

          # Get the exon number

          exonEndNumber = end_exon

          splitIndex = object$blockCount >= exonEndNumber

          if(TRUE %in% splitIndex){

            object = object[splitIndex]

            plusStrandRanges <- subset(object, strand(object) == "+")

            minusStrandRanges <- subset(object, strand(object) == "-")

          } else {

          "There are no genes with that many exons. Please enter another exon end."

          }
        } else {


          stop("If end_feature is set to Exon then end_exon must not be null")

        }
      }

    } else {

      stop("If start feature is TES then no exon start information must be entered.")

    }

  }


#----AT THIS POINT WE HAVE A PLUS STRAND GRANGES OBJECT, AND A MINUS STRAND GRANGES OBJECT. FILTERED FOR GENES--------------
  # THAT HAVE THE REQUIRED NUMBER OF EXONS. THEY MIGHT NOT HAVE ANY GENES.

  # Set the tempRanges start and end points to 0

  minusStrandTempRanges = minusStrandRanges
  plusStrandTempRanges = plusStrandRanges

  # If start feature is TSS, set the minus strand start points to the TSS (end).

  if(start_feature == 'TSS'){

    # Set start flanking region for strands

    if(start_direction == 'up'){

      start(plusStrandTempRanges) = start(plusStrandRanges) - start_flank

      start(minusStrandTempRanges) = 0

      end(minusStrandTempRanges) = end(minusStrandRanges) + start_flank

    } else {

      start(plusStrandTempRanges) = start(plusStrandRanges) + start_flank

      start(minusStrandTempRanges) = 0

      end(minusStrandTempRanges) = end(minusStrandRanges) - start_flank

    }

    # Set start flanking regions for minus strand

  }

  # If start feature is TES, set the plus strand start points to the TES (end)

  if(start_feature == 'TES'){

    if(start_direction == 'up'){

      start(plusStrandTempRanges) = end(plusStrandRanges) - start_flank

      end(minusStrandTempRanges) = start(minusStrandRanges) + start_flank

    } else {

      start(plusStrandTempRanges) = end(plusStrandRanges) + start_flank

      end(minusStrandTempRanges) = start(minusStrandRanges) + start_flank

    }

  }

  if(start_feature == 'Exon'){

    if(start_exon_boundary == '5prime'){

      #Set the plus_strand exon starts

      tss_plus <- start(plusStrandRanges)

      block_counts_plus <- plusStrandRanges$blockCount

      block_starts_plus <- lapply(strsplit(plusStrandRanges$blockStarts, ","), as.numeric)

      exon_starts_plus <- vector(mode = 'numeric', length = length(plusStrandRanges))

      for (i in seq_along(plusStrandRanges)) {
        exon_starts_plus[i] = tss_plus[i] + block_starts_plus[[i]][start_exon]
      }

      #start(plusStrandTempRanges) = exon_starts_plus

      #Set the minus_strand exon starts

      start_minus = start(minusStrandRanges)

      block_counts_minus = minusStrandRanges$blockCount

      block_starts_minus = lapply(strsplit(minusStrandRanges$blockStarts, ","), as.numeric)
      block_starts_minus = revElements(block_starts_minus,)
      block_sizes_minus = lapply(strsplit(minusStrandRanges$blockSizes, ","), as.numeric)
      block_sizes_minus = revElements(block_sizes_minus,)

      exon_starts_minus = vector(mode = 'numeric', length = length(minusStrandRanges))

      for (i in seq_along(minusStrandRanges)) {
        exon_starts_minus[i] <- start_minus[i] + block_starts_minus[[i]][start_exon] + block_sizes_minus[[i]][start_exon]
      }

      #start(minusStrandTempRanges) = 0
      #end(minusStrandTempRanges) = exon_starts_minus

    } else {

      # Get the 3' exon start regions

      tss_plus <- start(plusStrandRanges)
      block_counts_plus <- plusStrandRanges$blockCount

      block_sizes_plus = lapply(strsplit(plusStrandRanges$blockSizes, ","), as.numeric)
      block_starts_plus <- lapply(strsplit(plusStrandRanges$blockStarts, ","), as.numeric)

      exon_starts_plus <- vector(mode = 'numeric', length = length(plusStrandRanges))

      for (i in seq_along(plusStrandRanges)) {
        exon_starts_plus[i] <- tss_plus[i] + block_starts_plus[[i]][start_exon] + block_sizes_plus[[i]][start_exon]
      }

      #start(plusStrandTempRanges) = exon_starts_plus

      #Set the minus_strand exon starts

      start_minus = start(minusStrandRanges)

      block_counts_minus = minusStrandRanges$blockCount

      block_starts_minus = lapply(strsplit(minusStrandRanges$blockStarts, ","), as.numeric)
      block_starts_minus = revElements(block_starts_minus,)
      block_sizes_minus = lapply(strsplit(minusStrandRanges$blockSizes, ","), as.numeric)
      block_sizes_minus = revElements(block_sizes_minus,)

      exon_starts_minus = vector(mode = 'numeric', length = length(minusStrandRanges))

      for (i in seq_along(minusStrandRanges)) {
        exon_starts_minus[i] <- start_minus[i] + block_starts_minus[[i]][start_exon]
      }

      #start(minusStrandTempRanges) = 0
      #end(minusStrandTempRanges) = exon_starts_minus

    }

    if(start_direction == 'up'){

      start(plusStrandTempRanges) = exon_starts_plus - start_flank

      start(minusStrandTempRanges) = 0

      end(minusStrandTempRanges) = exon_starts_minus + start_flank

    } else {

      start(plusStrandTempRanges) = exon_starts_plus + start_flank

      start(minusStrandTempRanges) = 0

      end(minusStrandTempRanges) = exon_starts_minus - start_flank


    }

  }

#-----------AT THIS POINT WE HAVE ALL THE START POINTS SET CORRECTLY IN THE PLUS & MINUS TEMP RANGES OBJECTS--------------


  if(end_feature == 'TSS'){

    if(end_direction == 'up'){

      end(plusStrandTempRanges) = start(plusStrandRanges) - end_flank

      start(minusStrandTempRanges) = end(minusStrandRanges) + end_flank


    } else {

      end(plusStrandTempRanges) = start(plusStrandRanges) + end_flank

      start(minusStrandTempRanges) = end(minusStrandRanges) - end_flank

    }

  }


  if(end_feature == 'TES'){

    if(end_direction == 'up'){

      end(plusStrandTempRanges) = end(plusStrandRanges) - end_flank

      start(minusStrandTempRanges) = start(minusStrandRanges) + end_flank


    } else {

      end(plusStrandTempRanges) = end(plusStrandRanges) + end_flank

      start(minusStrandTempRanges) = start(minusStrandRanges) - end_flank

    }


  }

  else if (end_feature=="Exon"){

    if(end_exon_boundary == '5prime'){

      #Set the plus_strand exon ends

      tss_plus <- start(plusStrandRanges)

      block_counts_plus <- plusStrandRanges$blockCount

      block_starts_plus <- lapply(strsplit(plusStrandRanges$blockStarts, ","), as.numeric)

      exon_starts_plus <- vector(mode = 'numeric', length = length(plusStrandRanges))

      for (i in seq_along(plusStrandRanges)) {
        exon_starts_plus[i] = tss_plus[i] + block_starts_plus[[i]][end_exon]
      }

      #end(plusStrandTempRanges) = exon_starts_plus

      #Set the minus_strand exon ends

      start_minus = start(minusStrandRanges)

      block_counts_minus = minusStrandRanges$blockCount

      block_starts_minus = lapply(strsplit(minusStrandRanges$blockStarts, ","), as.numeric)
      block_starts_minus = revElements(block_starts_minus,)
      block_sizes_minus = lapply(strsplit(minusStrandRanges$blockSizes, ","), as.numeric)
      block_sizes_minus = revElements(block_sizes_minus,)

      exon_starts_minus = vector(mode = 'numeric', length = length(minusStrandRanges))

      for (i in seq_along(minusStrandRanges)) {
        exon_starts_minus[i] <- start_minus[i] + block_starts_minus[[i]][end_exon] + block_sizes_minus[[i]][end_exon]
      }

      #start(minusStrandTempRanges) = exon_starts_minus

    } else {

      # Get the 3' exon start regions

      tss_plus <- start(plusStrandRanges)
      block_counts_plus <- plusStrandRanges$blockCount

      block_sizes_plus = lapply(strsplit(plusStrandRanges$blockSizes, ","), as.numeric)
      block_starts_plus <- lapply(strsplit(plusStrandRanges$blockStarts, ","), as.numeric)

      exon_starts_plus <- vector(mode = 'numeric', length = length(plusStrandRanges))

      for (i in seq_along(plusStrandRanges)) {
        exon_starts_plus[i] <- tss_plus[i] + block_starts_plus[[i]][end_exon] + block_sizes_plus[[i]][end_exon]
      }

      #end(plusStrandTempRanges) = exon_starts_plus

      #Set the minus_strand exon starts

      start_minus = start(minusStrandRanges)

      block_counts_minus = minusStrandRanges$blockCount

      block_starts_minus = lapply(strsplit(minusStrandRanges$blockStarts, ","), as.numeric)
      block_starts_minus = revElements(block_starts_minus,)
      block_sizes_minus = lapply(strsplit(minusStrandRanges$blockSizes, ","), as.numeric)
      block_sizes_minus = revElements(block_sizes_minus,)

      exon_starts_minus = vector(mode = 'numeric', length = length(minusStrandRanges))

      for (i in seq_along(minusStrandRanges)) {
        exon_starts_minus[i] <- start_minus[i] + block_starts_minus[[i]][end_exon]
      }

      #start(minusStrandTempRanges) = exon_starts_minus

    }

    if(end_direction == 'up'){

      end(plusStrandTempRanges) = exon_starts_plus - end_flank

      start(minusStrandTempRanges) = exon_starts_minus + end_flank


    } else {

      end(plusStrandTempRanges) = exon_starts_plus + end_flank

      start(minusStrandTempRanges) = exon_starts_minus - end_flank

    }

  }

  if(return == "all"){

    allRanges = c(minusStrandTempRanges, plusStrandTempRanges)

    return(allRanges)

  }

  else if(return == "plus"){

    return(plusStrandTempRanges)

  }

  else if(return == "minus"){

    return(minusStrandTempRanges)

  }

  else{

    print("Please set a return value to 'all', 'plus' or 'minus'")

  }

}
