## functions to create, update and manipulate the sqlite database for meta-analysis 

#' @importFrom utils read.table
validateProjDataFile <- function(projDataFile, conn) {
  projData <- read.table(file = projDataFile, header = TRUE, 
                         sep = '\t', stringsAsFactors = FALSE)
  if(sum(grepl(pattern = '-', x = projData$sampleName)) > 0){
    warning("Converting dashes (-) in sample names to underscore (_)\n")
    projData$sampleName <- gsub(pattern = '-', '_', x = projData$sampleName)
  }
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
  
  #check 4 (when an existing database is updated, samples in projData object that clash with the
  #existing sample names will be ignored)
  if(RSQLite::dbExistsTable(conn, 'processedSamples')) {
    ps <- RSQLite::dbReadTable(conn, 'processedSamples')
    existingSamples <- intersect(ps$sampleName, projData$sampleName)
    if(length(existingSamples) > 0) {
      warning("The following samples already exists in the database. 
              These samples will be ignored during the database update:",
              paste0(existingSamples, collapse = ','))
    }
    projData <- rbind(ps, projData[!projData$sampleName %in% existingSamples,])
  }
  return(projData)
}

insertTableGtfData <- function(conn, name, gtfFilePath) {
  if(RSQLite::dbExistsTable(conn, name)) {
    message("Reading GTF data from the database")
    gtfData <- GenomicRanges::makeGRangesFromDataFrame(
      df = RSQLite::dbReadTable(conn, name), keep.extra.columns = TRUE)
  } else {
    message("Importing GTF annotations")
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
}

insertTableBedData <- function(conn, projData) {
  if(RSQLite::dbExistsTable(conn = conn, name = 'bedData')) {
    message("Reading existing BED data from database")
    bedData <- GenomicRanges::makeGRangesListFromDataFrame(
      df = RSQLite::dbReadTable(conn, 'bedData'), 
      split.field = 'sampleName')
    
    #find which samples in projData don't have any data yet
    newSamples <- setdiff(projData$sampleName, names(bedData))
    if(length(newSamples) > 0) {
      message("Reading BED files for new samples")
      newBedData <- importBedFiles(
        filePaths = projData[projData$sampleName %in% newSamples,]$bedFilePath,
        colnames = c('chrom', 'start', 'end', 'strand'))
      names(newBedData) <- newSamples
      
      message("Appending new interval datasets in 'bedData' table")
      dt <- as.data.table(newBedData)
      colnames(dt)[2] <- 'sampleName'
      RSQLite::dbWriteTable(conn = conn, name = 'bedData', 
                            value = dt, 
                            append = TRUE)
      bedData <- c(bedData, newBedData)
    } 
    return(bedData)
    } else {
    bedData <- importBedFiles(
      filePaths = projData$bedFilePath,
      colnames = c('chrom', 'start', 'end', 'strand'))
    names(bedData) <- projData$sampleName
    message("Saving interval datasets in 'bedData' table")
    dt <- as.data.table(bedData)
    colnames(dt)[2] <- 'sampleName'
    RSQLite::dbWriteTable(conn = conn, name = 'bedData', 
                          value = dt, 
                          append = TRUE)
    return(bedData)
  } 
}

insertTableAnnotationSummaries <- function(conn, bedData, ...) {
  
  if(RSQLite::dbExistsTable(conn, 'annotationSummaries')) {
    annotationSummaries <- RSQLite::dbReadTable(conn, 'annotationSummaries')
    
    newSamples <- setdiff(names(bedData), annotationSummaries$sampleName)
    if(length(newSamples) > 0) {
      message("Calculating annotation summaries for new samples")
      newAnnot <- summarizeQueryRegionsMulti(queryRegionsList = bedData[newSamples], 
                                             ...)
      newAnnot <- data.table::as.data.table(newAnnot, keep.rownames = TRUE)
      colnames(newAnnot)[1] <- 'sampleName'
      
      message("Appending new annotation summaries in 'annotationSummaries' table")
      RSQLite::dbWriteTable(conn = conn, name = 'annotationSummaries', 
                            value = newAnnot, 
                            append = TRUE)
    } 
  } else {
    message("Calculating annotation summaries")
    #get feature overlap table
    annotationSummaries <- summarizeQueryRegionsMulti(queryRegionsList = bedData, 
                                                      ...)
    annotationSummaries <- data.table::as.data.table(x = annotationSummaries, 
                                                     keep.rownames = TRUE)
    colnames(annotationSummaries)[1] <- 'sampleName'
    message("Saving annotation summaries in 'annotationSummaries' table")
    RSQLite::dbWriteTable(conn = conn, name = 'annotationSummaries', 
                          value = annotationSummaries, 
                          append = TRUE)
    
  }
}

insertTableFeatureBoundaryCoverageProfiles <- function(conn, bedData, ...) {
  
  if(RSQLite::dbExistsTable(conn, 'featureBoundaryCoverageProfiles')) {
    cvg <- RSQLite::dbReadTable(conn, 'featureBoundaryCoverageProfiles')
    
    newSamples <- setdiff(names(bedData), cvg$sampleName)
    if(length(newSamples) > 0) {
      message("Calculating feature boundary coverage profiles for new samples")
      newCvg <- getFeatureBoundaryCoverageMulti(bedData[newSamples], ...)

      message("Appending new coverage profiles in 
              'featureBoundaryCoverageProfiles' table")
      RSQLite::dbWriteTable(conn = conn, 
                            name = 'featureBoundaryCoverageProfiles', 
                            value = newCvg, 
                            append = TRUE)
    } 
  } else {
    message("Calculating feature boundary coverage profiles")
    cvg <- getFeatureBoundaryCoverageMulti(bedData, ...)
    message("Saving coverage profiles in 
              'featureBoundaryCoverageProfiles' table")
    RSQLite::dbWriteTable(conn = conn, 
                          name = 'featureBoundaryCoverageProfiles', 
                          value = cvg, 
                          append = TRUE)
    
  }
}

insertTableDiscoveredMotifs <- function(conn, bedData, txdbFeatures, ...) {
  
  if(RSQLite::dbExistsTable(conn, 'discoveredMotifs')) {
    motifData <- RSQLite::dbReadTable(conn, 'discoveredMotifs')
    
    newSamples <- setdiff(names(bedData), motifData$sampleName)
    if(length(newSamples) > 0) {
      message("Started motif discovery module for new samples...")
      newMotifData <- lapply(X = names(bedData[newSamples]),
                             FUN = function(s) {
                               queryRegions <- bedData[[s]]
                               message("Discovering motifs for sample ",s)
                               results <- discoverFeatureSpecificMotifs(
                                 queryRegions = queryRegions, 
                                 txdbFeatures = txdbFeatures, 
                                 ...)
                               results$sampleName <- s
                               return(results)
                             })
      newMotifData <- do.call(rbind, newMotifData)
      
      message("Appending new motif discovery stats to 
              'discoveredMotifs' table")
      RSQLite::dbWriteTable(conn = conn, 
                            name = 'discoveredMotifs', 
                            value = newMotifData, 
                            append = TRUE)
    } 
  } else {
    message("Started motif discovery module...")
    motifData <- lapply(X = names(bedData),
                        FUN = function(s) {
                          queryRegions <- bedData[[s]]
                          message("Discovering motifs for sample ",s)
                          results <- discoverFeatureSpecificMotifs(
                            queryRegions = queryRegions, 
                            txdbFeatures = txdbFeatures, 
                            ...)
                          results$sampleName <- s
                          return(results)
                        })
    motifData <- do.call(rbind, motifData)
    message("Saving motif discovery stats in   
            'discoveredMotifs' table")
    RSQLite::dbWriteTable(conn = conn, 
                          name = 'discoveredMotifs', 
                          value = motifData, 
                          append = TRUE)
  }
}

insertTableOverlapMatrix <- function(conn, name, bedData, ...) {
  if(RSQLite::dbExistsTable(conn, name)) {
    M <- RSQLite::dbReadTable(conn, name)
    
    newSamples <- setdiff(names(bedData), colnames(M))
    if(length(newSamples) > 0) {
      newM <- getIntervalOverlapMatrix(queryRegionsList = bedData[newSamples], ...)
      M <- cbind(data.table::as.data.table(M), 
                 data.table::as.data.table(newM))
    }
    RSQLite::dbWriteTable(conn = conn, name = name, 
                          value = M, 
                          overwrite = TRUE)
  } else {
      message("Running function: getIntervalOverlapMatrix for table",name)
      M <- getIntervalOverlapMatrix(queryRegionsList = bedData, ...)
      RSQLite::dbWriteTable(conn = conn, name = name, 
                            value = data.table::as.data.table(M, 
                                                              keep.rownames = TRUE), 
                            overwrite = TRUE)
  }
}

#' summarizeDatabaseContent 
#' 
#' Given a path to an sqlite database created using RCAS::createDB function,
#' accesses the database and provides a quick summary of available samples 
#' and number of entries of each sample in the available tables of the 
#' database.  
#' 
#' @param dbPath Path to the sqlite database
#' @return A data.frame object 
#' @import RSQLite
#' @export
summarizeDatabaseContent <- function(dbPath) {
  mydb <- RSQLite::dbConnect(RSQLite::SQLite(), dbPath)
  tables <- RSQLite::dbListTables(mydb)
  #exclude gtfData because it is independent from sample names
  tables <- tables[!tables %in% c('gtfData')]

  summary <- lapply(tables, function(tbl) {
    df <- RSQLite::dbReadTable(mydb, tbl)
    #tables that don't have 'sampleName' as a column have to be 
    #handled separately. 
    if (tbl == 'geneOverlaps') {
      r <- data.frame(apply(df[-1], 2, sum))
      colnames(r) <- tbl
      r$sampleName <- rownames(r)
      return(r)
    } else {
      r <- data.frame(table(df$sampleName))
      colnames(r) <- c('sampleName', tbl)
      return(r)
    }
  })
  summary <- Reduce(function(...) merge(..., all=T, by = 'sampleName'), summary)
  RSQLite::dbDisconnect(mydb)
  return(summary)
}

#' deleteSampleDataFromDB
#' 
#' Given a list of sample names, the function deletes all datasets calculated
#' for the given samples from the database.
#' 
#' @param dbPath Path to the sqlite database
#' @param sampleNames The names of the samples for which all relevant datasets
#'   should be deleted from the database. Tip: Use RSQLite::dbReadTable function
#'   to read the table 'processedSamples' to see which samples are available in
#'   the database.
#' @return SQLiteConnection object with updated contents in the dbPath 
#' @import RSQLite
#' @export
deleteSampleDataFromDB <- function(dbPath, sampleNames) {
  mydb <- RSQLite::dbConnect(RSQLite::SQLite(), dbPath)
  tables <- RSQLite::dbListTables(mydb)
  #exclude gtfData because it is independent from sample names
  tables <- tables[!tables %in% c('gtfData')]
  
  for(tbl in tables) {
    #tables that don't have 'sampleName' as a column have to be 
    #handled separately. tables with the column 'sampleName' 
    #can be updated all at the same time (see the else statement)
    if(tbl == 'geneOverlaps') {
      df <- RSQLite::dbReadTable(mydb, tbl)
      df <- df[,!colnames(df) %in% sampleNames]
      RSQLite::dbRemoveTable(mydb, 'geneOverlaps')
      RSQLite::dbWriteTable(mydb, 'geneOverlaps', df, append = TRUE)
      dfNew <- RSQLite::dbReadTable(mydb, 'geneOverlaps')
      message("Removing columns '",paste0(sampleNames, 
                                          collapse = ', '),
              "' from table ",tbl,
              "... Columns affected: ",ncol(dfNew)-ncol(df))
    } else {
      rs <- RSQLite::dbSendStatement(conn = mydb, 
                                     statement = paste0(
                                       'DELETE FROM ',tbl,
                                       ' WHERE "sampleName" = :s'))
      RSQLite::dbBind(rs, params = list(s = sampleNames))
      message("Removing rows '",paste0(sampleNames, 
                                       collapse = ', '),
              "' from table ",tbl,
              "... Rows affected: ",RSQLite::dbGetRowsAffected(rs))
      dbClearResult(rs)
    }
  }
  return(mydb)
}


#' createDB
#' 
#' Creates an sqlite database consisting of various tables of data obtained 
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
#' @param update TRUE/FALSE (default: FALSE) whether an existing database 
#'   should be updated 
#' @param genomeVersion  A character string to denote for which genome version 
#'   the analysis is being done. Available options are hg19/hg38 (human), mm9/mm10 
#'   (mouse), ce10 (worm) and dm3 (fly).
#' @param annotationSummary TRUE/FALSE (default:TRUE) whether annotation summary
#' module should be run
#' @param coverageProfiles TRUE/FALSE (default: TRUE) whether coverage profiles
#' module should be run
#' @param motifAnalysis TRUE/FALSE (default: TRUE) whether motif discovery module
#' should be run 
#' @param nodeN Number of cpus to use for parallel processing (default: 1)
#' @return Path to an SQLiteConnection object created by RSQLite package
#' @importFrom RSQLite dbConnect
#' @importFrom RSQLite dbWriteTable
#' @importFrom RSQLite dbDisconnect
#' @importFrom RSQLite SQLite
#' @importFrom data.table as.data.table
#' @importFrom proxy dist
#' @examples 
#' FUS_path <- system.file("extdata", "FUS_Nakaya2013c_hg19.bed",
#' package='RCAS') 
#' 
#' FMR1_path <- system.file("extdata",
#' "FMR1_Ascano2012a_hg19.bed", package='RCAS') 
#' 
#' projData <- data.frame('sampleName' = c('FUS', 'FMR1'), 
#' 'bedFilePath' = c(FUS_path,FMR1_path), stringsAsFactors = FALSE)
#'  
#' write.table(projData, 'myProjDataFile.tsv', sep = '\t', quote =FALSE, 
#' row.names = FALSE) 
#' 
#' gtfFilePath <- system.file("extdata",
#' "hg19.sample.gtf", package='RCAS')
#' 
#' createDB(dbPath = 'hg19.RCASDB.sqlite', 
#' projDataFile = './myProjDataFile.tsv', 
#' gtfFilePath = gtfFilePath,  
#' genomeVersion = 'hg19', 
#' motifAnalysis = FALSE, 
#' coverageProfiles = FALSE)
#' 
#' #Note: to add new data to an existing database, set update = TRUE
#' 
#' @export
createDB <- function(dbPath = file.path(getwd(), 'rcasDB.sqlite'), 
                     projDataFile, gtfFilePath = '', update = FALSE, 
                     genomeVersion, 
                     annotationSummary = TRUE, 
                     coverageProfiles = TRUE, 
                     motifAnalysis = TRUE, 
                     nodeN = 1) {
  
  #check genome version
  if (!genomeVersion %in% c('hg19', 'hg38', 'dm3', 'ce10', 'mm9', 'mm10')){
    stop (genomeVersion,' is not a supported genome version.')
  }
  
  if(file.exists(dbPath) & update == FALSE) {
    stop("A database already exists at ",dbPath,"\n",
         "Either provide a different path or set argument 'update' to TRUE  
         in order to update this database")
  } else if(update == FALSE & !file.exists(gtfFilePath)) {
    stop("A GTF file doesn't exist at:\n",gtfFilePath,"\n")
  } else if(!file.exists(projDataFile)) {
    stop("A project data file doesn't exist at:",projDataFile,"\n")
  }
  #create connection to sqlite 
  mydb <- RSQLite::dbConnect(RSQLite::SQLite(), dbPath)
  #read and validate project data file 
  projData <- validateProjDataFile(projDataFile, mydb) 
  #import and save genome annotations
  gtfData <- insertTableGtfData(conn = mydb, name = 'gtfData',
                                gtfFilePath = gtfFilePath)
  message("Parsing transcript features")
  txdbFeatures <- getTxdbFeaturesFromGRanges(gtfData)
  
  bedData <- insertTableBedData(conn = mydb,projData = projData) 

  if(annotationSummary == TRUE) {
    insertTableAnnotationSummaries(conn = mydb, bedData = bedData, 
                                   txdbFeatures = txdbFeatures, nodeN = nodeN)
    
    insertTableOverlapMatrix(conn = mydb, name = 'geneOverlaps',
                             bedData = bedData, 
                             targetRegions = gtfData[gtfData$type == 'gene',],
                             targetRegionNames = gtfData[gtfData$type == 'gene',]$gene_name,
                             nodeN = nodeN)
  }
  
  if(coverageProfiles == TRUE) {
    insertTableFeatureBoundaryCoverageProfiles(conn = mydb, bedData = bedData, 
                                             txdbFeatures = txdbFeatures, 
                                             sampleN = 10000)
  }
  
  if(motifAnalysis == TRUE) {
    insertTableDiscoveredMotifs(conn = mydb, 
                                bedData = bedData, 
                                txdbFeatures = txdbFeatures, 
                                genomeVersion = genomeVersion, 
                                motifN = 1,
                                nCores = 2,
                                resizeN = 15, 
                                sampleN = 10000)
  }
  
  RSQLite::dbWriteTable(mydb, 'processedSamples', projData, overwrite = TRUE)
  RSQLite::dbDisconnect(mydb)
  message("Finished preparing the database and disconnected.",
      "\nType",paste0(" RSQLite::dbConnect(RSQLite::SQLite(), \'",
                    dbPath,"\')")," to reconnect")
  return(dbPath)
} 
