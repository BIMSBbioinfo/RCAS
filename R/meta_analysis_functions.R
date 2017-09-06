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

insertTableAnnotationSummaries <- function(conn, ...) {
  cat("Calculating annotation summaries\n")
  #get feature overlap table
  annotationSummaries <- summarizeQueryRegionsMulti(...)
  
  cat("Saving annotation summaries in 'annotationSummaries' table\n")
  RSQLite::dbWriteTable(conn = conn, name = 'annotationSummaries', 
                        value = as.data.frame(annotationSummaries), 
                        append = TRUE)
}

insertTableOverlapMatrix <- function(conn, name, ...) {
  cat("Running function: getIntervalOverlapMatrix for table",name,"\n")
  M <- getIntervalOverlapMatrix(...)

  RSQLite::dbWriteTable(conn = conn, name = name, 
                        value = data.table::as.data.table(M, 
                                                          keep.rownames = TRUE), 
                        append = TRUE)
}

insertTableBedData <- function(conn, sampleNames, ...) {
  cat("Reading bed files\n")
  bedData <- importBedFiles(...)
  names(bedData) <- sampleNames
  
  cat("Saving interval datasets in 'bedData' table\n")
  RSQLite::dbWriteTable(conn = conn, name = 'bedData', 
                        value = as.data.table(bedData), 
                        append = TRUE)
  return(bedData)
}

validateProjDataFile <- function(projDataFile) {
  projData <- read.table(file = projDataFile, header = TRUE, 
                         sep = '\t', stringsAsFactors = FALSE)
  
  #check 1
  ambigousSamples <- projData[table(projData$sampleName) > 1,]$sampleName
  if(length(ambigousSamples) > 0) {
    stop("Some sample names are ambiguously assigned to multiple files in:\n",projDataFile,"\n",
         "See sample names:",paste0(ambigousSamples, collapse = ', '),"\n")
  }
  #check 2
  ambigousFiles <- projData[table(projData$bedFilePath) > 1,]$bedFilePath
  if(length(ambigousFiles) > 0) {
    stop("Some bed files are ambiguously assigned to multiple sample names",projDataFile,"\n",
         "See:",paste0(ambigousFiles, collapse = ', '),"\n")
  }
  
  #check 3
  unaccessibleFiles <- projData$bedFilePath[!sapply(projData$bedFilePath, 
                                                    file.exists)]
  if(length(unaccessibleFiles) > 0) {
    stop("Following files are not accessible or do not exist at the provided location:\n",
         "See:",paste0(unaccessibleFiles, collapse = ', '),"\n")
  }
  return(projData)
}

#' createDB
#' 
#' Creates and sqlite database consisting of various tables of data obtained 
#' from processed BED files
#' 
#' @param dbPath Path to the sqlite database file (could be an existing file or 
#'   a new file path to be created at the given path)
#' @param projDataFile A file consisting of meta-data about the input samples. 
#'   Must minimally consist of two columns: 1. sampleName (name of the sample) 
#'   2. bedFilePath (full path to the location of the BED file containing data 
#'   for the sample)
#' @param gtfFilePath Path to the GTF file (preferably downloaded from the 
#'   Ensembl database) that contains genome annotations
#' @param nodeN Number of cpus to use for parallel processing (default: 1)
#' @return Path to an SQLiteConnection object created by RSQLite package
#' @importFrom RSQLite dbConnect
#' @importFrom RSQLite dbWriteTable
#' @importFrom RSQLite dbDisconnect
#' @importFrom RSQLite SQLite
#' @importFrom data.table as.data.table
#' @importFrom proxy dist
#' @export
createDB <- function(dbPath = file.path(getwd(), 'rcasDB.sqlite'), 
                     projDataFile, gtfFilePath, nodeN = 1) {
  if(file.exists(dbPath)) {
    stop("A database already exists at ",dbPath,"\n",
         "Either provide a different path or use updateDB() 
         function to update an existing database")
  } else if(!file.exists(gtfFilePath)) {
    stop("A GTF file doesn't exist at:\n",gtfFilePath,"\n")
  } else if(!file.exists(projDataFile)) {
    stop("A project data file doesn't exist at:",projDataFile,"\n")
  }
  #create connection to sqlite 
  mydb <- RSQLite::dbConnect(RSQLite::SQLite(), dbPath)
  #read and validate project data file
  projData <- validateProjDataFile(projDataFile)
  cat("Importing GTF annotations\n")
  gffData <- importGtf(filePath = gtfFilePath)
  cat("Parsing transcript features\n")
  txdbFeatures <- getTxdbFeaturesFromGRanges(gffData)
  
  bedData <- insertTableBedData(conn = mydb, 
                                sampleNames = projData$sampleName,
                                filePaths = projData$bedFilePath,
                                colnames = c('chrom', 'start', 'end', 'strand'))
  insertTableAnnotationSummaries(conn = mydb, queryRegionsList = bedData,
                          txdbFeatures = txdbFeatures, nodeN = nodeN)
  insertTableOverlapMatrix(conn = mydb, name = 'geneOverlaps',
                           queryRegionsList = bedData, 
                           targetRegions = gffData[gffData$type == 'gene',],
                           targetRegionNames = gffData[gffData$type == 'gene',]$gene_name,
                           nodeN = nodeN)
  RSQLite::dbWriteTable(mydb, 'processedSamples', projData)
  RSQLite::dbDisconnect(mydb)
  cat("Finished preparing the database and disconnected.\n",
      "Type",paste0("RSQLite::dbConnect(RSQLite::SQLite(), \'",
                    dbPath,"\')"),"to reconnect\n")
  return(dbPath)
} 











