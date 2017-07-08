#' importBedFiles
#'
#' This function is a wrapper that uses \code{RCAS::importBed()} function to 
#' import BED files as a GRangesList object
#' 
#' @param filePaths A vector of paths to one or more BED files
#' @param ... Other parameters passed to RCAS::importBed and
#'   rtracklayer::import.bed function
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
  bedData <- GenomicRanges::GRangesList(sapply(filePaths, function(f) {
    RCAS::importBed(filePath = f, ...)
  }, USE.NAMES = TRUE))
  names(bedData) <- gsub(pattern = '.bed$', 
                         replacement = '', 
                         x = basename(names(bedData)), 
                         ignore.case = TRUE)
  return(bedData)
}


#' summarizeQueryRegionsMulti
#' 
#' This function is a wrapper function to run RCAS::summarizeQueryRegions 
#' multiple times, which is useful to get a matrix of overlap counts between a 
#' list of BED files with a txdbFeatures extracted from GTF file
#' 
#' @param queryRegionsList GRangesList object imported from multiple BED files
#'   using \code{importBedFiles} function
#' @param txdbFeatures List of GRanges objects - outputs of 
#'   \code{getTxdbFeaturesFromGRanges} and \code{getTxdbFeatures} functions
#' @return A list consisting of two data.frame objects: one for raw overlap 
#'   counts and one for percentage of overlap counts (raw overlap counts divided
#'   by the number of query regions in the corresponding BED file)
#' @examples
#' data(gff)
#' data(queryRegions)
#' queryRegionsList <- GRangesList(queryRegions, queryRegions)
#' names(queryRegionsList) <- c('q1', 'q2')
#' txdbFeatures <- getTxdbFeaturesFromGRanges(gffData = gff)
#' summaryMatrix <- summarizeQueryRegionsMulti(queryRegionsList = queryRegionsList,
#'                                  txdbFeatures = txdbFeatures)
#'
#' @import GenomicRanges
#' @importFrom S4Vectors queryHits
#' @export
summarizeQueryRegionsMulti <- function(queryRegionsList, txdbFeatures) {
  summaryRaw <- sapply(queryRegionsList, function(x) summarizeQueryRegions(x, txdbFeatures))
  rownames(summaryRaw) <- names(txdbFeatures)
  return(t(summaryRaw))
}

#' getOverlapSimilarityMatrix
#' 
#' This function takes as input a GRangesList objects and finds pairwise
#' overlapping intervals in all list elements. Then, it creates a similarity
#' matrix based on the percentage of intervals in each GRanges object to overlap
#' any intervals in every other GRanges object in the list.
#' 
#' @param queryRegionsList A GRangesList object
#' @param ... Other arguments passed to IRanges::overlapsAny function
#' @return A matrix with equal number of rows and columns where the size is 
#'   determined by the number of list elements in the input object. Matrix
#'   values are between 0 and 1, where 0 means no similarity (no overlap) and 1
#'   means all elements of object A overlap have at least one overlap with some
#'   element in B (but not necessarily vice versa).
#' @examples 
#' data(queryRegions)
#' queryRegionsList <- GRangesList(queryRegions, queryRegions, 
#'                                 queryRegions, queryRegions)
#' getOverlapSimilarityMatrix(queryRegionsList, type = 'any')                                
#' @importFrom IRanges overlapsAny
#' @export
getOverlapSimilarityMatrix <- function(queryRegionsList, ...) {
  sapply(queryRegionsList, function(gr1) {
    sapply(queryRegionsList, function(gr2) {
      round(sum(IRanges::overlapsAny(gr1, gr2, ...)) /
              length(gr1), 2)
    })
  })
}

getOverlapSimilarityMatrix2 <- function(queryRegionsList, ...) {
  sapply(queryRegionsList, function(gr1) {
    sapply(queryRegionsList, function(gr2) {
      round(sum(IRanges::overlapsAny(gr1, gr2, ...)) /
              length(gr1), 2)
    })
  })
}


