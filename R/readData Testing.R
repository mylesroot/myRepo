myBedFile = ("~/Documents/Shaun Data/Saccharomyces_cerevisiae.EF4.74_SGDv64_CUTandSUT_withUTRs_noEstimates_antisense_intergenic_4xlncRNAs_final.pyCheckGTFfile.bed")

myBed = readBed(myBedFile)

myGtfFile = ("~/Documents/Shaun Data/Saccharomyces_cerevisiae.EF4.74_SGDv64_CUTandSUT_withUTRs_noEstimates_antisense_intergenic_4xlncRNAs_final.pyCheckGTFfile.output.gtf")

myGtf = readGFFAsGRanges(myGtfFile)

proteinCodingGeneList = subset(myGtf, myGtf$"gene_biotype" == "protein_coding")$gene_id


#--------------------------------------------------------Filter By Introns------------------------------------------------------------


#extract introns and put them into a GRanges object

filterForIntrons = function(bed){
  bed = subset(bed,bed$blockCount>1)
  return(bed)

}



#--------------------------------------------------------Filter By Gene List------------------------------------------------------------
rpGenes = "~/Documents/Shaun Data/RP_Intron_Genes.bed"
rpGenes = readBed(rpGenes)
rpGeneNames = rpGenes$name

# This list can be used to specify genes that are misannotated or to subset the GRanges by

genesOfInterest = c("tK(UUU)P", "RDN37-1", "YLR316C", "YPL198W", "Q0045", "YNL162W-A", "YDL130W-A")

# Check the a list of genes and select for/select inverse from the GRanges object

# inverse parameter to get the inverse of the filter

filter_genes <- function(bed, geneList, inverse = FALSE) {
  # Check if genes in geneList are present in the GRanges object
  found_genes <- bed$name[bed$name %in% geneList]
  not_found_genes <- geneList[!geneList %in% bed$name]

  if (length(not_found_genes) > 0) {
    message("The following genes were not found in the GRanges object: ",
            paste(not_found_genes, collapse = ", "))
  }

  if (length(found_genes) > 0) {
    message("Genes from the list that were found in the GRanges object: ", paste(found_genes, collapse = ", "))


    if(inverse){
      proceed <- readline(prompt = "Proceed to remove found genes? (y/n): ")
    }

    else {
      proceed <- readline(prompt = "Proceed to remove all genes but the one/s found? (y/n): ")
    }

    if (tolower(proceed) == "y") {

      if (inverse) {
        # Filter to REMOVE the genes in the list
        filtered_bed <- bed[!bed$name %in% geneList, ]
        message("These genes were removed ", paste(found_genes, collapse = ", "))
      } else {
        # Filter to KEEP the specified genes
        filtered_bed <- bed[bed$name %in% geneList, ]
        message("GRanges object filtered to the specified genes.")
      }
    } else {
      message("Filtering cancelled.")
      filtered_bed = bed
    }

  } else {

    message("No genes from the list were found, filtering cancelled")
    filtered_bed = bed

  }

  return(filtered_bed)

}

#--------------------------------------------------------Filter By Strand------------------------------------------------------------


filter_by_strand <- function(bed, plus = FALSE, minus = FALSE, other = FALSE) {

  # Check if at least one strand type is selected for filtering
  if (!plus && !minus && !other) {
    stop("At least one strand type must be selected for filtering.")
  }

  keep_strands <- c()
  if (plus) keep_strands <- c(keep_strands, "+")
  if (minus) keep_strands <- c(keep_strands, "-")
  if (other) keep_strands <- c(keep_strands, "*")

  filtered_bed <- subset(bed, strand(bed) %in% keep_strands)

  return(filtered_bed)
}

#--------------------------------------------------------Get Ranges Function---------------------------------------------------------

getRanges = function(object, plus = FALSE, minus = FALSE, upstream = 2000, downstream = 2000, exon1 = TRUE, intron1 = TRUE, exon2Tes = TRUE){

  # Ensure input validity
  if (!inherits(object, "GRanges")) {
    stop("Input must be a GRanges object.")
  }

  # Create GRanges objects to concatenate
  upstreamRange = GRanges()
  downstreamRange = GRanges()
  exon1Range = GRanges()
  intron1Range = GRanges()
  exon2TesRange = GRanges()

  keep_strands <- c()
  if (plus) keep_strands <- c(keep_strands, "+", "*")
  if (minus) keep_strands <- c(keep_strands, "-")

  object <- subset(object, strand(object) %in% keep_strands)

  # If upstream != NULL, get upstream flank regions and append to list

  if(!is.null(upstream)){

    upstreamRange = promoters(object, upstream = upstream, downstream = 0)

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
    tss_exon1_ranges <- GRanges(
      seqnames = c(seqnames(plus_strand), seqnames(minus_strand)),
      ranges = IRanges(c(start(plus_strand), exon1_ends_minus),
                       c(exon1_ends_plus, end(minus_strand))),
      strand = c(rep("+", length(plus_strand)), rep("-", length(minus_strand)))
    )

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
    intron1_ranges <- GRanges(
      seqnames = c(seqnames(plus_strand), seqnames(minus_strand)),
      ranges = IRanges(c(intron1_starts_plus, intron1_ends_minus),
                       c(intron1_ends_plus, intron1_starts_minus)),
      strand = c(rep("+", length(plus_strand)), rep("-", length(minus_strand)))
    )

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
    tes_exon2start_ranges <- GRanges(
      seqnames = c(seqnames(plus_strand), seqnames(minus_strand)),
      ranges = IRanges(c(exon2_starts_plus, exon2_starts_minus),
                       c(end(plus_strand), tss_minus)),
      strand = c(rep("+", length(plus_strand)), rep("-", length(minus_strand)))
    )

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
    TES_downstream <- GRanges(
      seqnames = c(seqnames(plus_strand), seqnames(minus_strand)),
      ranges = IRanges(c(tes_plus, downstream_minus),
                       c(downstream_plus, tes_minus)),
      strand = c(rep("+", length(plus_strand)), rep("-", length(minus_strand)))
    )

    # Add metadata (optional)
    if (length(mcols(object)) > 0) {
      mcols(TES_downstream) <- mcols(object)
    }

    downstreamRange = TES_downstream

  }


  # Create a GRanges list of all the ranges. If the range is empty, then remove it from the final list.

  ranges = GRangesList(upstreamRange, exon1Range, intron1Range, exon2TesRange, downstreamRange)
  keep_elements = list()

  for(i in seq_along(ranges)){

    keep_elements = c(keep_elements, length(ranges[[i]])>0)

  }

  ranges = ranges[keep_elements==TRUE]

  return(ranges)

}
