#' Filter GRanges Objects With Block Counts >1
#'
#' This function takes a GRanges object and then returns a new GRanges object
#' where all the ranges have more than 1 block count. I.e. they contain introns.
#'
#' @param object The GRanges object that you would like to filter
#'
#' @export


filterForIntrons = function(object){
  object = subset(object,object$blockCount>1)
  return(object)

}


#' Filter GRanges Objects Using A Predefined List Of Gene Names
#'
#' This function takes a character vector of gene names and filters the bed file
#' to either exclude, or include the genes from the list.
#'
#' @param object The GRanges object that you would like to filter
#' @param geneList A character vector of gene names
#' @param inverse A TRUE/FALSE variable to specify whether you want to return the inverse of the list provided
#'
#' @export


filterByList <- function(object, geneList, inverse = FALSE) {
  # Check if genes in geneList are present in the GRanges object
  found_genes <- object$name[object$name %in% geneList]
  not_found_genes <- geneList[!geneList %in% object$name]

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
        filtered_object <- object[!object$name %in% geneList, ]
        message("These genes were removed ", paste(found_genes, collapse = ", "))
      } else {
        # Filter to KEEP the specified genes
        filtered_object <- object[object$name %in% geneList, ]
        message("GRanges object filtered to the specified genes.")
      }
    } else {
      message("Filtering cancelled.")
      filtered_object = object
    }

  } else {

    message("No genes from the list were found, filtering cancelled")
    filtered_object = object

  }

  return(filtered_object)

}

#--------------------------------------------------------Filter By Strand------------------------------------------------------------

#' Filter a GRanges object by strandedness of the ranges
#'
#' Returns a GRanges object with only the ranges of the strandedness you want
#'
#' @param object The GRanges object that you would like to filter
#' @param plus TRUE/FALSE whether you want to return the plus strand ranges
#' @param minus TRUE/FALSE whether you want to return the minus strand ranges
#' @param other TRUE/FALSE whether you want to return ranges with no strandedness - marked '*'
#'
#' @export

filterByStrand <- function(object, plus = FALSE, minus = FALSE, other = FALSE) {

  # Check if at least one strand type is selected for filtering
  if (!plus && !minus && !other) {
    stop("At least one strand type must be selected for filtering.")
  }

  keep_strands <- c()
  if (plus) keep_strands <- c(keep_strands, "+")
  if (minus) keep_strands <- c(keep_strands, "-")
  if (other) keep_strands <- c(keep_strands, "*")

  filtered_bed <- subset(object, strand(object) %in% keep_strands)

  return(filtered_bed)
}
