#' importBedFiles
#'
#' This function is a wrapper that uses \code{RCAS::importBed()} function to
#' import BED files as a GRangesList object
#' 
#' @param filePaths A vector of paths to one or more BED files
#' @param ... Parameters passed to RCAS::importBed function
#' @return A \code{GRangesList} object containing the coordinates of the
#'   intervals from multiple input BED files
#'   
#' @examples
#' input1 <- system.file("extdata", "testfile.bed", package='RCAS')
#' input2 <- system.file("extdata", "testfile2.bed", package='RCAS')
#' importBedFiles(filePaths = c(input1, input2), keepStandardChr = TRUE)
#' 
#' @importFrom GenomicRanges GRangesList
#' @export
importBedFiles <- function(filePaths, ...) {
  GenomicRanges::GRangesList(sapply(filePaths, function(f) {
    RCAS::importBed(filePath = f, ...)
  }, USE.NAMES = TRUE))
}

