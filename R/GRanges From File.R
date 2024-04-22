#' Load A High .bed or a .gtf file
#'
#' This function takes the file path of a .bed or .gtf file and returns a
#' GRanges object of the information in the file
#'
#' @param filePath A file path to the .bed or .gtf file
#' @param gtf TRUE/FALSE if the file is a .gtf
#' @param bed TRUE/FALSE if the file is a .bed
#' @importFrom rtracklayer readGFFAsGRanges
#' @importFrom genomation readBed
#' @export
#'

grangesFromFile = function(filePath, gtf = FALSE, bed = FALSE){

  if(gtf){

    file = rtracklayer::readGFFAsGRanges(filePath)

  }

  if(bed){

    file = genomation::readBed(filePath)

  }

  return(file)

}
