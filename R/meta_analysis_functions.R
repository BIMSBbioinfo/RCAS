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
#' bedData <- importBedFiles(filePaths = c(input1, input2), 
#' keepStandardChr = TRUE) 
#' # when importing multiple bed files with different column names, it 
#' # is required to pass the common column names to be parsed from the 
#' # bed files
#' bedData <- importBedFiles(filePaths = c(input1, input2),
#'                  colnames = c('chrom', 'start', 'end', 'strand'))
#' 
#' @importFrom GenomicRanges GRangesList
#' @importFrom pbapply pbsapply
#' @export
importBedFiles <- function(filePaths, ...) {
  bedData <- GenomicRanges::GRangesList(pbapply::pbsapply(filePaths, function(f) {
    RCAS::importBed(filePath = f, debug = FALSE, ...)
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
#' @param nodeN Positive integer value that denotes the number of cpus to use
#'   for parallel processing (default: 1)
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
#' @importFrom parallel makeCluster
#' @importFrom parallel clusterExport
#' @importFrom parallel stopCluster
#' @importFrom pbapply pbsapply 
#' @export
summarizeQueryRegionsMulti <- function(queryRegionsList, txdbFeatures, nodeN = 1) {
  cl <- parallel::makeCluster(nodeN)
  parallel::clusterExport(cl = cl, 
                          varlist = c('txdbFeatures', 'summarizeQueryRegions'), 
                          envir = environment())
  summaryRaw <- pbapply::pbsapply(queryRegionsList, 
                                  function(x) {
                                    summarizeQueryRegions(x, txdbFeatures)
                                    }, cl = cl)
  rownames(summaryRaw) <- c(names(txdbFeatures), 'NoFeatures')
  parallel::stopCluster(cl)
  return(t(summaryRaw))
}

#' getFeatureBoundaryCoverageMulti
#' 
#' This function is a wrapper function to run RCAS::getFeatureBoundaryCoverage 
#' multiple times, which is useful to get coverage signals across different 
#' kinds of transcript features for a given list of bed files imported as a 
#' GRangesList object.
#' 
#' @param bedData GRangesList object imported from multiple BED files 
#'   using \code{importBedFiles} function
#' @param txdbFeatures List of GRanges objects - outputs of 
#'   \code{getTxdbFeaturesFromGRanges} and \code{getTxdbFeatures} functions
#' @param sampleN (default=10000) Positive integer value that is used to 
#'   randomly down-sample the target feature coordinates to improve the runtime.
#'   Set to 0 to avoid downsampling.
#'   
#' @return A data.frame object with coverage data at three prime and five prime
#'   boundaries of a list of transcript features
#' @examples
#' data(gff)
#' data(queryRegions)
#' queryRegionsList <- GRangesList(queryRegions, queryRegions)
#' names(queryRegionsList) <- c('q1', 'q2')
#' txdbFeatures <- getTxdbFeaturesFromGRanges(gffData = gff)
#' getFeatureBoundaryCoverageMulti(queryRegionsList, txdbFeatures, sampleN = 500)
#' 
#' @importFrom pbapply pblapply
#' @export
getFeatureBoundaryCoverageMulti <- function(bedData, 
                                            txdbFeatures, 
                                            sampleN = 10000) {
  results <- pblapply(X = names(bedData),
    FUN = function(s) {
      queryCoords <- bedData[[s]]
      cvg <- lapply(names(txdbFeatures),
                    function(featureType) {
                      featureCoords <- txdbFeatures[[featureType]]
                      cvgF <- getFeatureBoundaryCoverage(
                        queryRegions = queryCoords,
                        featureCoords = featureCoords,
                        sampleN = sampleN,
                        boundaryType = 'fiveprime')
                      
                      cvgT <- getFeatureBoundaryCoverage(
                        queryRegions = queryCoords,
                        featureCoords = featureCoords,
                        sampleN = sampleN,
                        boundaryType = 'threeprime')
                      
                      cvgF$boundary <- 'fiveprime'
                      cvgT$boundary <- 'threeprime'
                      
                      df <- rbind(cvgF, cvgT)
                      df$feature <- featureType
                      return(df)
                      })
      cvg <- do.call(rbind, cvg)
      cvg$sampleName <- s
      return(cvg)
      })
  return(do.call(rbind, results))
}


#' getIntervalOverlapMatrix
#' 
#' This function is used to obtain a binary matrix of overlaps between a list of
#' GRanges objects (GRangesList object) and a target GRanges object. The
#' resulting matrix has N rows where N is the number of intervals in the target
#' GRanges object and M columns where M is the number GRanges objects in the
#' query GRangesList object.
#' 
#' @param queryRegionsList A GRangesList object
#' @param targetRegions A GRanges object
#' @param targetRegionNames Optional vector of names to be used as rownames in
#'   the resulting matrix. The vector indices must correspond to the intervals
#'   in targetRegions object.
#' @param nodeN Positive integer value to use one or more cpus for parallel
#'   computation (default: 1)
#' @return A binary matrix object consisting of number of rows equal to the
#'   number of intervals in targetRegions object, and number of columns equal to
#'   the number of GRanges objects available in the queryRegionsList object.
#' @examples 
#' data(gff)
#' input1 <- system.file("extdata", "testfile.bed", package='RCAS') 
#' input2 <- system.file("extdata", "testfile2.bed", package='RCAS') 
#' bedData <- RCAS::importBedFiles(filePaths = c(input1, input2)) 
#' M <- RCAS::getIntervalOverlapMatrix(
#' queryRegionsList = bedData, 
#' targetRegions = gff[gff$type == 'gene',][1:100], 
#' targetRegionNames = gff[gff$type == 'gene',][1:100]$gene_name)
#' @export
getIntervalOverlapMatrix <- function(queryRegionsList, targetRegions, targetRegionNames = NULL, nodeN = 1) {
  
  if(!is.null(targetRegionNames) & length(targetRegions) != length(targetRegionNames)) {
    stop("The sizes of targetRegions and targetRegionNames must be equal\n")
  }
  cl <- parallel::makeCluster(nodeN)
  parallel::clusterExport(cl = cl, 
                          varlist = c('targetRegions'), 
                          envir = environment())
  summaryRaw <- pbapply::pbsapply(queryRegionsList, 
                                  function(x) {
                                    myOverlaps <- GenomicRanges::findOverlaps(targetRegions, x)
                                    unique(S4Vectors::queryHits(myOverlaps))
                                  }, cl = cl)
  parallel::stopCluster(cl)
  
  rowN <- length(targetRegions)
  colN <- length(queryRegionsList)
  
  M <- sapply(X = summaryRaw, 
              FUN = function(x) {
                v <- rep(0, rowN) 
                v[x] <- 1
                return(v)
              })
  if(!is.null(targetRegionNames)) {
    rownames(M) <- targetRegionNames
  }
  return(M)
}


#' runReportMetaAnalysis
#' 
#' Generate a stand-alone HTML report with interactive figures and tables from a
#' pre-calculated RCAS database (using RCAS::createDB) to compare multiple 
#' samples.
#' 
#' @param dbPath Path to the sqlite database generated by RCAS::createDB
#' @param sampleTablePath A tab-separated file with two columns (no rownames) 
#'   header 1: sampleName, header 2: sampleGroup
#' @param outDir Path to the output directory. (default: current working 
#'   directory)
#' @param outFile Name of the output HTML report (by default, the base name of
#'   sampleTablePath value is used to create a name for the HTML report)
#' @param quiet boolean value (default: FALSE). If set to TRUE, progress bars 
#'   and chunk labels will be suppressed while knitting the Rmd file.
#' @param selfContained boolean value (default: TRUE). By default, the generated
#'   html file will be self-contained, which means that all figures and tables 
#'   will be embedded in a single html file with no external dependencies (See 
#'   rmarkdown::html_document)
#' @return An html generated using rmarkdown/knitr/pandoc that contains 
#'   interactive figures, tables, and text that provide an overview of the 
#'   experiment
#' @examples
#' dbPath <- system.file("extdata", "hg19.RCASDB.sqlite",
#' package='RCAS') 
#' 
#' #Hint: use RCAS::summarizeDatabaseContent to see which samples have processed
#' #data in the database. 
#' summarizeDatabaseContent(dbPath = dbPath)
#' 
#' #Create a data table for samples and their groups sampleGroup field is used
#' #to group replicates of #the same sample into one group in visualizations. 
#' #Any arbitrary name can be used for sampleGroup field. However, entries in
#' #the sampleName field must be available in the queried database
#' sampleData <- data.frame('sampleName' = c('FUS', 'FMR1'), 
#'     'sampleGroup' = c('FUS', 'FMR1'), stringsAsFactors = FALSE) 
#' write.table(sampleData, 'sampleDataTable.tsv', sep = '\t', 
#'     quote =FALSE, row.names = FALSE)
#' #Use the generated database to run a report 
#' runReportMetaAnalysis(dbPath = 'hg19.RCASDB.sqlite', 
#'     sampleTablePath = 'sampleDataTable.tsv')
#' @import rmarkdown
#' @import knitr
#' @importFrom plotly ggplotly
#' @import ggplot2
#' @importFrom pheatmap pheatmap
#' @importFrom cowplot plot_grid
#' @importFrom ggseqlogo geom_logo
#' @importFrom ggseqlogo theme_logo
#' @import DT
#' @export
runReportMetaAnalysis <- function(dbPath = 'RCAS.sqlite',
                                  sampleTablePath,
                                  outDir = getwd(),
                                  outFile = NULL,
                                  quiet = FALSE,
                                  selfContained = TRUE) {
  
  reportFile <- system.file("reporting_scripts", 
                            "reportMetaAnalysis.Rmd", package='RCAS')
  headerFile <- system.file("reporting_scripts", 
                            "header.html", package='RCAS')
  footerFile <- system.file("reporting_scripts", 
                            "footer.html", package='RCAS')
  if(is.null(outFile)) {
    outFile <- paste0(basename(sampleTablePath), '.RCAS.metaAnalysisReport.html')
  } 
  rmarkdown::render(
    input = reportFile, 
    output_dir = outDir,
    intermediates_dir = outDir,
    output_file = outFile,
    output_format = rmarkdown::html_document(
      toc = TRUE,
      toc_float = TRUE,
      theme = 'simplex',
      number_sections = TRUE,
      includes = rmarkdown::includes(in_header = headerFile, 
                                     after_body = footerFile), 
      self_contained = selfContained
    ),
    params = list(dbPath = dbPath,
                  sampleTablePath = sampleTablePath,
                  workdir = outDir),
    quiet = quiet
  )
}
