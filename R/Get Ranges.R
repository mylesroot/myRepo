#' Get Specific Gene Ranges From Your GRanges Object
#'
#' This function returns a list of GRanges objects corresponding to specific gene ranges. Included in that list
#' is the upstream range, TSS->Exon1 range, Intron1 range, Exon2->TES range and the downstream range. This GRanges list
#' can then be used by the regionPlot and plotRegion functions to plot these specific regions on a graph. IMPORTANT: You must
#' filter your GRanges object to only include intron-containing genes before using this function
#'
#' @param object The GRanges object that you would like to get ranges from
#' @param plus TRUE/FALSE to only return ranges with plus strandedness
#' @param minus TRUE/FALSE to only return ranges with minus strandedness
#' @param upstream The width of the upstream range you want to return. Enter NULL if you don't want the upstream range.
#' @param downstream The width of the downstream range you want to return. Enter NULL if you don't want the downstream range.
#' @param exon1 TRUE/FALSE to return the TSS->Exon1 range in your list
#' @param intron1 TRUE/FALSE to return the intron1 range in your list
#' @param exon2Tes TRUE/FALSE to return the exon2 -> TES range
#'
#' @importFrom GenomicRanges GRanges
#' @importFrom GenomicRanges GRangesList
#' @importFrom IRanges promoters
#' @import BiocGenerics
#'
#'@export


getRanges = function(object, plus = TRUE, minus = FALSE, upstream = 2000, downstream = 2000, exon1 = TRUE, intron1 = TRUE, exon2Tes = TRUE){

  # Ensure input validity
  if (!inherits(object, "GRanges")) {
    stop("Input must be a GRanges object.")
  }

  # Create GRanges objects to concatenate
  upstreamRange = GenomicRanges::GRanges()
  downstreamRange = GenomicRanges::GRanges()
  exon1Range = GenomicRanges::GRanges()
  intron1Range = GenomicRanges::GRanges()
  exon2TesRange = GenomicRanges::GRanges()

  keep_strands <- c()
  if (plus) keep_strands <- c(keep_strands, "+", "*")
  if (minus) keep_strands <- c(keep_strands, "-")

  object <- subset(object, BiocGenerics::strand(object) %in% keep_strands)

  # If upstream != NULL, get upstream flank regions and append to list

  if(!is.null(upstream)){

    upstreamRange = IRanges::promoters(object, upstream = upstream, downstream = 0)

  }

  # If exon1 is set to TRUE, get exon1 regions and append to list

  if(exon1){

    # Separate strands for processing
    plus_strand <- subset(object, strand(object) == "+")
    minus_strand <- subset(object, strand(object) == "-")

    # Calculate TSS-Exon1 for plus strand
    tss_plus <- start(plus_strand)
    exon_counts_plus <- plus_strand$blockCount

    exon_sizes_plus <- lapply(strsplit(plus_strand$blockSizes, ","), as.numeric)
    exon_starts_plus <- lapply(strsplit(plus_strand$blockStarts, ","), as.numeric)

    exon1_ends_plus <- vector(mode = 'numeric', length = length(plus_strand))
    for (i in seq_along(plus_strand)) {
      exon1_ends_plus[i] <- tss_plus[i] + exon_starts_plus[[i]][1] + exon_sizes_plus[[i]][1]
    }

    # Calculate TSS-Exon1 for minus strand (working backwards)
    tss_minus <- end(minus_strand)  # End becomes the effective TSS
    exon_counts_minus <- minus_strand$blockCount
    tes_minus = start(minus_strand)

    # COULD ADD: Error check that block count = no. of block sizes

    exon_sizes_minus <- (lapply(strsplit(minus_strand$blockSizes, ","), as.numeric))
    exon_sizes_minus = revElements(exon_sizes_minus,)# reverse order
    exon_starts_minus <- (lapply(strsplit(minus_strand$blockStarts, ","), as.numeric))
    exon_starts_minus = revElements(exon_starts_minus,)# reverse order

    exon1_ends_minus <- vector(mode = 'numeric', length = length(minus_strand))
    for (i in seq_along(minus_strand)) {
      exon1_ends_minus[i] <- tes_minus[i] + exon_starts_minus[[i]][1]
    }

    # Combine results and create new GRanges object

    tss_exon1_ranges = GenomicRanges::GRanges(seqnames = c(seqnames(plus_strand), seqnames(minus_strand)),
                                              ranges = IRanges(c(start(plus_strand), exon1_ends_minus),
                                                               c(exon1_ends_plus, end(minus_strand))),
                                              strand = c(rep("+", length(plus_strand)), rep("-", length(minus_strand))))

    # Add metadata (optional)
    if (length(mcols(object)) > 0) {
      mcols(tss_exon1_ranges) <- mcols(object)
    }

    exon1Range = tss_exon1_ranges

  }


  # If intron1 is set to TRUE, get intron1 regions and append to list

  if(intron1){

    # Separate strands for processing
    plus_strand <- subset(object, strand(object) == "+")
    minus_strand <- subset(object, strand(object) == "-")

    # Calculate Exon1 ends (intron 1 starts) for plus strand
    tss_plus <- start(plus_strand)
    exon_counts_plus <- plus_strand$blockCount

    exon_sizes_plus <- lapply(strsplit(plus_strand$blockSizes, ","), as.numeric)
    exon_starts_plus <- lapply(strsplit(plus_strand$blockStarts, ","), as.numeric)

    intron1_starts_plus <- vector(mode = 'numeric', length = length(plus_strand))
    for (i in seq_along(plus_strand)) {
      intron1_starts_plus[i] <- tss_plus[i] + exon_starts_plus[[i]][1] + exon_sizes_plus[[i]][1]
    }


    # Calculate exon1 ends (intron 1 starts) for minus strand (working backwards)
    tss_minus <- end(minus_strand)  # End becomes the effective TSS
    exon_counts_minus <- minus_strand$blockCount
    tes_minus = start(minus_strand)

    # COULD ADD: Error check that block count = no. of block sizes

    exon_sizes_minus <- (lapply(strsplit(minus_strand$blockSizes, ","), as.numeric))
    exon_sizes_minus = revElements(exon_sizes_minus,)# reverse order
    exon_starts_minus <- (lapply(strsplit(minus_strand$blockStarts, ","), as.numeric))
    exon_starts_minus = revElements(exon_starts_minus,)# reverse order

    intron1_starts_minus <- vector(mode = 'numeric', length = length(minus_strand))
    for (i in seq_along(minus_strand)) {
      intron1_starts_minus[i] <- tes_minus[i] + exon_starts_minus[[i]][1]
    }

    # Calculate exon 2 starts (intron 1 ends) for plus strand

    intron1_ends_plus <- vector(mode = 'numeric', length = length(plus_strand))
    for (i in seq_along(plus_strand)) {
      intron1_ends_plus[i] <- tss_plus[i] + exon_starts_plus[[i]][2]
    }


    # Calculate exon 2 starts (intron 1 ends) for minus strand

    intron1_ends_minus <- vector(mode = 'numeric', length = length(minus_strand))
    for (i in seq_along(minus_strand)) {
      intron1_ends_minus[i] <- tes_minus[i] + exon_starts_minus[[i]][2] + exon_sizes_minus[[i]][2]
    }

    # Combine results and create new GRanges object

    intron1_ranges = GenomicRanges::GRanges(seqnames = c(seqnames(plus_strand), seqnames(minus_strand)),
                                            ranges = IRanges(c(intron1_starts_plus, intron1_ends_minus),
                                                             c(intron1_ends_plus, intron1_starts_minus)),
                                            strand = c(rep("+", length(plus_strand)), rep("-", length(minus_strand))))

    # Add metadata (optional)
    if (length(mcols(object)) > 0) {
      mcols(intron1_ranges) <- mcols(object)
    }

    intron1Range = intron1_ranges

  }


  # If exon2Tes is set to TRUE, get exon2Tes regions and append to list

  if(exon2Tes) {

    # Separate strands for processing
    plus_strand <- subset(object, strand(object) == "+")
    minus_strand <- subset(object, strand(object) == "-")

    # Calculate Exon2 starts for plus strand

    tss_plus <- start(plus_strand)
    exon_counts_plus <- plus_strand$blockCount

    exon_sizes_plus <- lapply(strsplit(plus_strand$blockSizes, ","), as.numeric)
    exon_starts_plus <- lapply(strsplit(plus_strand$blockStarts, ","), as.numeric)

    exon2_starts_plus <- vector(mode = 'numeric', length = length(plus_strand))
    for (i in seq_along(plus_strand)) {
      exon2_starts_plus[i] <- tss_plus[i] + exon_starts_plus[[i]][2]
    }

    # Calculate exon 2 starts (intron 1 ends) for minus strand

    tss_minus <- end(minus_strand)  # End becomes the effective TSS
    exon_counts_minus <- minus_strand$blockCount
    tes_minus = start(minus_strand)

    exon_sizes_minus <- (lapply(strsplit(minus_strand$blockSizes, ","), as.numeric))
    exon_sizes_minus = revElements(exon_sizes_minus,)# reverse order
    exon_starts_minus <- (lapply(strsplit(minus_strand$blockStarts, ","), as.numeric))
    exon_starts_minus = revElements(exon_starts_minus,)# reverse order

    exon2_starts_minus <- vector(mode = 'numeric', length = length(minus_strand))
    for (i in seq_along(minus_strand)) {
      exon2_starts_minus[i] <- tes_minus[i] + exon_starts_minus[[i]][2] + exon_sizes_minus[[i]][2]
    }

    # Combine results and create new GRanges object

    tes_exon2start_ranges = GenomicRanges::GRanges(seqnames = c(seqnames(plus_strand), seqnames(minus_strand)),
                                                   ranges = IRanges(c(exon2_starts_plus, exon2_starts_minus),
                                                                    c(end(plus_strand), tss_minus)),
                                                   strand = c(rep("+", length(plus_strand)), rep("-", length(minus_strand))))

    # Add metadata (optional)
    if (length(mcols(object)) > 0) {
      mcols(tes_exon2start_ranges) <- mcols(object)
    }

    exon2TesRange = tes_exon2start_ranges


  }

  # If downstream !=NULL, get downstream regions and append to list

  if(!is.null(downstream)){

    # Separate strands for processing
    plus_strand <- subset(object, strand(object) == "+")
    minus_strand <- subset(object, strand(object) == "-")

    tes_minus = start(minus_strand)
    tes_plus = end(plus_strand)

    downstream_plus = tes_plus + downstream
    downstream_minus = tes_minus - downstream

    # Combine results and create new GRanges object

    TES_downstream = GenomicRanges::GRanges(seqnames = c(seqnames(plus_strand), seqnames(minus_strand)),
                                            ranges = IRanges(c(tes_plus, downstream_minus),
                                                             c(downstream_plus, tes_minus)),
                                            strand = c(rep("+", length(plus_strand)), rep("-", length(minus_strand))))

    # Add metadata (optional)
    if (length(mcols(object)) > 0) {
      mcols(TES_downstream) <- mcols(object)
    }

    downstreamRange = TES_downstream

  }


  # Create a GRanges list of all the ranges. If the range is empty, then remove it from the final list.

  ranges = GenomicRanges::GRangesList(upstreamRange, exon1Range, intron1Range, exon2TesRange, downstreamRange)
  keep_elements = list()

  for(i in seq_along(ranges)){

    keep_elements = c(keep_elements, length(ranges[[i]])>0)

  }

  ranges = ranges[keep_elements==TRUE]

  return(ranges)

}
