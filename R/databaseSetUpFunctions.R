## functions to create, update and manipulate the sqlite database for meta-analysis 
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

insertTableGtfData <- function(conn, name, gtfFilePath) {
  cat("Importing GTF annotations\n")
  gtfData <- importGtf(filePath = gtfFilePath, 
                       readFromRds = FALSE, 
                       colnames = c('source', 'type', 'score', 
                                    'gene_id', 'gene_name', 
                                    'transcript_id', 'exon_id', 
                                    'exon_number'))
  RSQLite::dbWriteTable(conn = conn, name = name, 
                        value = data.table::as.data.table(gtfData))
  return(gtfData)
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
  #import and save genome annotations
  gtfData <- insertTableGtfData(conn = mydb, name = 'gtfData', 
                                gtfFilePath = gtfFilePath)
  cat("Parsing transcript features\n")
  txdbFeatures <- getTxdbFeaturesFromGRanges(gtfData)
  
  bedData <- insertTableBedData(conn = mydb, 
                                sampleNames = projData$sampleName,
                                filePaths = projData$bedFilePath,
                                colnames = c('chrom', 'start', 'end', 'strand'))
  insertTableAnnotationSummaries(conn = mydb, queryRegionsList = bedData,
                                 txdbFeatures = txdbFeatures, nodeN = nodeN)
  insertTableOverlapMatrix(conn = mydb, name = 'geneOverlaps',
                           queryRegionsList = bedData, 
                           targetRegions = gtfData[gtfData$type == 'gene',],
                           targetRegionNames = gtfData[gtfData$type == 'gene',]$gene_name,
                           nodeN = nodeN)
  RSQLite::dbWriteTable(mydb, 'processedSamples', projData)
  RSQLite::dbDisconnect(mydb)
  cat("Finished preparing the database and disconnected.\n",
      "Type",paste0("RSQLite::dbConnect(RSQLite::SQLite(), \'",
                    dbPath,"\')"),"to reconnect\n")
  return(dbPath)
} 
