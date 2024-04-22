calculate_intron1_range <- function(bed) {

  # ADD ERROR CHECKING: ARE THERE POSITIVE AND NEGATIVE STRANDS PRESENT?

  # Ensure input validity
  if (!inherits(bed, "GRanges")) {
    stop("Input must be a GRanges object.")
  }

  # Separate strands for processing
  plus_strand <- subset(bed, strand(bed) == "+")
  minus_strand <- subset(bed, strand(bed) == "-")
  print(plus_strand)
  print(minus_strand)

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
  print(exon_starts_minus)
  print(exon_sizes_minus)

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

  intron1_ends_minus <- vector(mode = 'numeric', length = length(plus_strand))
  for (i in seq_along(minus_strand)) {
    intron1_ends_minus[i] <- tes_minus[i] + exon_starts_minus[[i]][2] + exon_sizes_minus[[i]][2]
  }

  # Combine results and create new GRanges object
  intron1_ranges <- GRanges(
    seqnames = c(seqnames(plus_strand), seqnames(minus_strand)),
    ranges = IRanges(c(intron1_starts_plus, intron1_ends_minus),
                     c(intron1_ends_plus, intron1_starts_minus)),
    strand = c(rep("+", length(plus_strand)), rep("-", length(minus_strand)))
  )

  # Add metadata (optional)
  if (length(mcols(bed)) > 0) {
    mcols(intron1_ranges) <- mcols(bed)
  }

  return(intron1_ranges)
}

