#' Install dependecies for genomation
install.packages( c("data.table","plyr","reshape2","ggplot2","gridBase","devtools"))
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install(c("GenomicRanges","rtracklayer","impute","Rsamtools"))

#' install genomation
library(devtools)
install_github("BIMSBbioinfo/genomation",build_vignettes=FALSE)

library(genomation)

# install GenomicRanges

BiocManager::install("GenomicRanges")
library(GenomicRanges)

# Could create a function here to store the filenames of the data for easy use

myBedFile = ("~/Documents/Shaun Data/Saccharomyces_cerevisiae.EF4.74_SGDv64_CUTandSUT_withUTRs_noEstimates_antisense_intergenic_4xlncRNAs_final.pyCheckGTFfile.bed")

myBedData = readGeneric(myBedFile, chr=1, start=2, end=3, strand=6, meta.cols=list(name=4, score=5, thickStart=7, thickEnd=8, itemRGB=9, blockCount=10, blockSize=11, blockStart=12))
myBedData
myBed = readBed(myBedFile)

myBed

#extract introns and put them into a GRanges object

bed12Introns = function(bed){

  bed = subset(bed,bed$blockCount>1)
  return(bed)

}


misannotatedGeneList <- c("YNL162W-A", "YDL130W-A", "Test")

# Check the list of misannotated genes and remove them from the GRanges object

# -v parameter to get the inverse of the filter

filter_misannotated_genes <- function(bed, misannotatedGenes) {
  # Check for misannotated genes
  found_genes <- bed$name[bed$name %in% misannotatedGenes]
  not_found_genes <- misannotatedGenes[!misannotatedGenes %in% bed$name]

  if (length(found_genes) > 0) {
    message("Misannotated genes present: ", paste(found_genes, collapse = ", "))

  if (length(not_found_genes) > 0) {
    message("Genes not found: ", paste(not_found_genes, collapse = ", "))
  }

    # Ask for user input
    proceed <- readline(prompt = "Proceed with filtering? (y/n): ")

    if (tolower(proceed) == "y") {
      # Filter the GRanges object
      filtered_bed <- bed[!bed$name %in% misannotatedGenes, ]
      return(filtered_bed)

    } else {
      message("Filtering cancelled.")
      return(bed) # Return the original bed object
    }

  } else {
    message("No misannotated genes found.")
    return(bed) # Return the original bed object
  }
}

# Filter for the plus strands only

plus_strand_filtered = function(bed){

  bed = subset(bed, strand(bed) =="+")

}



# Extracting 200kb upstream of the TSS from a GRanges object and creating a new GRanges the range and NO METADATA

upstream = function(gr, upstream=200, downstream=0, use.names=TRUE){

  gr = promoters(gr, upstream, downstream, use.names)

  mcols(gr) = NULL

  return(gr)

}

# Extracting the range of TSS->Start of Exon 1 and creating a new GRanges object for it with NO METADATA

# Take out any exon (Exon1 = TSS->Exon1 End)

calculate_tss_exon1_range <- function(gr) {

  # Ensure input validity
  if (!inherits(gr, "GRanges")) {
    stop("Input must be a GRanges object.")
  }

  # Extract necessary information
  tss <- start(gr)

  exon_counts <- gr$blockCount

  exon_sizes = lapply(strsplit(gr$blockSizes, ","), as.numeric)
  # exon_sizes <- strsplit(granges(granges_object)$blockSizes, ",")

  # Check block counts and block sizes match up

  exon_starts = lapply(strsplit(gr$blockStarts, ","), as.numeric)
  # exon_starts <- strsplit(granges(granges_object)$blockStarts, ",")

  # Calculate end of Exon 1
  # exon1_ends <- tss + exon_starts[[1]] + exon_sizes[[1]] - 1

  # Calculate end of Exon 1 (iterating over each feature)
  exon1_ends <- vector(mode = 'numeric', length = length(gr))
  for (i in seq_along(gr)) {
    exon1_ends[i] <- tss[i] + exon_starts[[i]][1] + exon_sizes[[i]][1]
  }

  # Create the new GRanges object
  tss_exon1_ranges <- GRanges(seqnames = seqnames(gr),
                              ranges = IRanges(start(gr), exon1_ends),
                              strand = '+'
                              )

  # Add metadata (optional)
  # if (length(mcols(granges_object)) > 0) {
  # mcols(tss_exon1_ranges) <- mcols(granges_object)
  # }

  return(tss_exon1_ranges)
}

# Create a new GRanges object with the start of intron 1 to the end of intron 1


calculate_intron1_range <- function(gr) {

  # Ensure input validity
  if (!inherits(gr, "GRanges")) {
    stop("Input must be a GRanges object.")
  }

  # Extract necessary information
  tss <- start(gr)

  exon_counts <- gr$blockCount

  exon_sizes = lapply(strsplit(gr$blockSizes, ","), as.numeric)
  # exon_sizes <- strsplit(granges(granges_object)$blockSizes, ",")

  exon_starts = lapply(strsplit(gr$blockStarts, ","), as.numeric)
  # exon_starts <- strsplit(granges(granges_object)$blockStarts, ",")

  # Calculate end of Exon 1
  # exon1_ends <- tss + exon_starts[[1]] + exon_sizes[[1]] - 1

  # Calculate end of Exon 1 (iterating over each feature)
  exon1_ends <- vector(mode = 'numeric', length = length(gr))
  for (i in seq_along(gr)) {
    exon1_ends[i] <- tss[i] + exon_starts[[i]][1] + exon_sizes[[i]][1] - 1
  }

  exon2_starts = vector(mode = 'numeric', length = length(gr))
  for (i in seq_along(gr)) {
    exon2_starts[i] = tss[i] + exon_starts[[i]][2] - 1
  }

  first_intron_ends = exon2_starts - 1

  first_intron_starts = exon1_ends + 1

  # Create the new GRanges object
  intron1_ranges <- GRanges(seqnames = seqnames(gr),
                              ranges = IRanges(first_intron_starts, first_intron_ends),
                              strand = '+'
                            )

  # Add metadata (optional)
  # if (length(mcols(granges_object)) > 0) {
  # mcols(tss_exon1_ranges) <- mcols(granges_object)
  # }

  return(intron1_ranges)
}



calculate_exon2_TES_range <- function(gr) {

  # Ensure input validity
  if (!inherits(gr, "GRanges")) {
    stop("Input must be a GRanges object.")
  }

  # Extract necessary information
  tes <- end(gr)

  exon_counts <- gr$blockCount

  exon_sizes = lapply(strsplit(gr$blockSizes, ","), as.numeric)
  # exon_sizes <- strsplit(granges(granges_object)$blockSizes, ",")

  exon_starts = lapply(strsplit(gr$blockStarts, ","), as.numeric)
  # exon_starts <- strsplit(granges(granges_object)$blockStarts, ",")

  # Calculate end of Exon 1
  # exon1_ends <- tss + exon_starts[[1]] + exon_sizes[[1]] - 1

  # Calculate the start of exon 2, iterating over each feature

  exon2_starts = vector(mode = 'numeric', length = length(gr))
  for (i in seq_along(gr)) {
    exon2_starts[i] = tss[i] + exon_starts[[i]][2] - 1
  }


  # Create the new GRanges object
  exon2_TES_ranges <- GRanges(seqnames = seqnames(gr),
                            ranges = IRanges(exon2_starts, tes),
                            strand = '+'
  )

  # Add metadata (optional)
  # if (length(mcols(granges_object)) > 0) {
  # mcols(tss_exon1_ranges) <- mcols(granges_object)
  # }

  return(exon2_TES_ranges)
}


# Calculating the TES to 200b downstream of TES

calculate_TES_200b_downstream = function(gr){

  tes = end(gr)

  downstream_200 = tes + 200

  # Create the new GRanges object
  TES_200b_downstream <- GRanges(seqnames = seqnames(gr),
                              ranges = IRanges(tes, downstream_200),
                              strand = '+'
  )

  return(TES_200b_downstream)

}


bedTest = upstream200(bedPlus)




