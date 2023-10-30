.datatable.aware=TRUE

#' importGtf
#'
#' This function uses \code{rtracklayer::import.gff()} function to import genome
#' annoatation data from an Ensembl gtf file
#'
#' @param filePath Path to a GTF file
#' @param saveObjectAsRds TRUE/FALSE (default:TRUE). If it is set to TRUE, a
#'   GRanges object will be created and saved in RDS format
#'   (<filePath>.granges.rds) so that importing can re-use this .rds file in
#'   next run.
#' @param readFromRds TRUE/FALSE (default:TRUE). If it is set to TRUE,
#'   annotation data will be imported from previously generated .rds file
#'   (<filePath>.granges.rds).
#' @param overwriteObjectAsRds TRUE/FALSE (default:FALSE). If it is set to TRUE,
#'   existing .rds file (<filePath>.granges.rds) will overwritten.
#' @param keepStandardChr TRUE/FALSE (default:TRUE). If it is set to TRUE,
#'   \code{seqlevelsStyle} will be converted to 'UCSC' and
#'   \code{keepStandardChromosomes} function  will be applied to only keep data
#'   from the standard chromosomes.
#' @param ... Other arguments passed to rtracklayer::import.gff function
#' @return A \code{GRanges} object containing the coordinates of the annotated
#'   genomic features in an input GTF file
#'
#' @examples
#' #import the data and write it into a .rds file
#' \dontrun{
#' importGtf(filePath='./Ensembl75.hg19.gtf')
#' }
#' #import the data but don't save it as RDS
#' \dontrun{
#' importGtf(filePath='./Ensembl75.hg19.gtf', saveObjectAsRds = FALSE)
#' }
#' #import the data and overwrite the previously generated
#' \dontrun{
#' importGtf(filePath='./Ensembl75.hg19.gtf', overwriteObjectAsRds = TRUE)
#' }
#'
#' @importFrom rtracklayer import.gff
#' @importFrom GenomeInfoDb seqlevelsStyle
#' @importFrom GenomeInfoDb keepStandardChromosomes
#' @export
importGtf <- function (filePath,  saveObjectAsRds = TRUE, readFromRds = TRUE,
                       overwriteObjectAsRds = FALSE, keepStandardChr = TRUE, ...) {

  rdsFilePath <- paste0(filePath, ".granges.rds")

  if (readFromRds == TRUE && file.exists(rdsFilePath)){
    message('Reading existing granges.rds object from ',rdsFilePath)
    gff <- readRDS(rdsFilePath)
  } else {
    message('importing gtf file from ', filePath)
    gff <- rtracklayer::import.gff(filePath, ...)
  }

  if (!is.null(gff)) {
    if (keepStandardChr == TRUE) {
      message('Keeping standard chromosomes only')
      GenomeInfoDb::seqlevelsStyle(gff) = 'UCSC'
      gff = GenomeInfoDb::keepStandardChromosomes(gff, pruning.mode = 'coarse')
    }

    if (saveObjectAsRds == TRUE) {
      if (!file.exists(rdsFilePath)){
        message('Saving gff object in RDS file at ',rdsFilePath)
        saveRDS(gff, file=rdsFilePath)
      } else {
        if (overwriteObjectAsRds == TRUE) {
          message('Overwriting gff object in RDS file at ',rdsFilePath)
          saveRDS(gff, file=rdsFilePath)
        } else {
          message('File ',rdsFilePath,' already exists.
                  Use overwriteObjectAsRds = TRUE to overwrite the file')
        }
      }
    return (gff)
    }
  } else{
    stop("Couldn't import the gtf from",filePath,"\n")
  }
}


#' importBed
#'
#' This function uses \code{rtracklayer::import.bed()} function to import BED
#' files
#'
#' @param filePath Path to a BED file
#' @param sampleN A positive integer value. The number of intervals in the input
#'   BED file are randomly downsampled to include intervals as many as
#'   \code{sampleN}. The input will be downsampled only if this value is larger
#'   than zero and less than the total number of input intervals.
#' @param keepStandardChr TRUE/FALSE (default:TRUE). If set to TRUE, will
#'   convert the \code{seqlevelsStyle} to 'UCSC' and apply
#'   \code{keepStandardChromosomes} function to only keep data from the standard
#'   chromosomes
#' @param debug TRUE/FALSE (default:TRUE). Set to FALSE to turn off messages
#' @param ... Other arguments passed to rtracklayer::import.bed function
#' @return A \code{GRanges} object containing the coordinates of the intervals
#'   from an input BED file
#'
#' @examples
#' input <- system.file("extdata", "testfile.bed", package='RCAS')
#' importBed(filePath = input, keepStandardChr = TRUE)
#'
#' @importFrom rtracklayer import.bed
#' @importFrom GenomeInfoDb seqlevelsStyle
#' @importFrom GenomeInfoDb keepStandardChromosomes
#' @export
importBed <- function (filePath, sampleN = 0, keepStandardChr = TRUE, debug = TRUE, ...) {

  if (file.exists(filePath)) {
    data = rtracklayer::import.bed(filePath, ...)
    if(debug == TRUE) {
      message('Processing ',filePath)
    }
    if (keepStandardChr == TRUE) {
      GenomeInfoDb::seqlevelsStyle(data) <- 'UCSC'
      if(debug == TRUE) {
        message('Keeping standard chromosomes only')
      }
      data <- GenomeInfoDb::keepStandardChromosomes(data, pruning.mode = 'coarse')
    }

    if (sampleN > 0 && sampleN < length(data)) {
      if(debug == TRUE) {
        message ('Downsampling intervals to size:',sampleN)
      }
      data <- data[sample(x = 1:length(data), size = sampleN)]
    }
    return(data)
  } else {
    stop("Cannot import BED file.
         It does not exist at the given location:",filePath,"\n")
  }
}


#' getTxdbFeatures
#' 
#' This function is deprecated. Use getTxdbFeaturesFromGRanges instead. 
#' 
#' @export
getTxdbFeatures <- function () {
  .Deprecated("getTxdbFeaturesFromGRanges")
}

#' getTxdbFeaturesFromGRanges
#'
#' This function takes as input a GRanges object that contains GTF file contents
#' (e.g from the output of \code{importGtf} function). 
#' Then extracts the coordinates of gene features
#' such as promoters, introns, exons, 5'/3' UTRs and
#' whole transcripts.
#'
#' @param gffData A GRanges object imported by \code{importGtf} function
#'
#' @examples
#' data(gff)
#' txdbFeatures <- getTxdbFeaturesFromGRanges(gffData = gff)
#'
#' @return A list of GRanges objects
#'
#' @import GenomicFeatures
#' @import GenomicRanges
#' @importFrom BiocGenerics unlist
#' @export
getTxdbFeaturesFromGRanges <- function (gffData) {

  txdb <- GenomicFeatures::makeTxDbFromGRanges(gffData)

  transcripts <- GenomicFeatures::transcripts(txdb)
  m <- match(transcripts$tx_name, gffData$transcript_id)
  transcripts$gene_name <- gffData[m]$gene_name

  tmp <- GenomicFeatures::exonsBy(x = txdb, by = "tx", use.names = TRUE)
  exons <- BiocGenerics::unlist(tmp)
  exons$tx_name <- names(exons)
  exons$gene_name <- gffData[match(names(exons), gffData$transcript_id)]$gene_name

  tmp <- GenomicFeatures::intronsByTranscript(x = txdb, use.names = TRUE)
  introns <- BiocGenerics::unlist(tmp)
  introns$tx_name <- names(introns)
  m <- match(names(introns), gffData$transcript_id)
  introns$gene_name <- gffData[m]$gene_name

  promoters <- GenomicFeatures::promoters(txdb)
  m <- match(promoters$tx_name, gffData$transcript_id)
  promoters$gene_name <- gffData[m]$gene_name

  #fiveUTRsByTranscript will return exonic regions of UTRs. These exonic regions
  #must be turned into a single region containing both exons and introns using
  #the range function we can get min-start and max-end of all UTR regions per
  #transcript
  tmp <- range(GenomicFeatures::fiveUTRsByTranscript(x = txdb, use.names = TRUE))
  fiveUTRs <- BiocGenerics::unlist(tmp)
  fiveUTRs$tx_name <- names(fiveUTRs)
  m <- match(names(fiveUTRs), gffData$transcript_id)
  fiveUTRs$gene_name <- gffData[m]$gene_name

  #threeUTRsByTranscript will return exonic regions of UTRs. These exonic regions
  #must be turned into a single region containing both exons and introns using
  #the range function we can get min-start and max-end of all UTR regions per
  #transcript
  tmp <- range(GenomicFeatures::threeUTRsByTranscript(x = txdb, use.names = TRUE))
  threeUTRs <- BiocGenerics::unlist(tmp)
  threeUTRs$tx_name <- names(threeUTRs)
  m <- match(names(threeUTRs), gffData$transcript_id)
  threeUTRs$gene_name <- gffData[m]$gene_name

  tmp <- GenomicFeatures::cdsBy(x = txdb, by = "tx", use.names = TRUE)
  cds <- BiocGenerics::unlist(tmp)
  cds$tx_name <- names(cds)
  cds$gene_name <- gffData[match(names(cds), gffData$transcript_id)]$gene_name

  txdbFeatures = list(
    'transcripts' = transcripts,
    'exons'       = exons,
    'promoters'   = promoters,
    'fiveUTRs'    = fiveUTRs,
    'introns'     = introns,
    'cds'         = cds,
    'threeUTRs'   = threeUTRs
  )
  return(txdbFeatures)
}

#' @importFrom data.table data.table
#' @importFrom S4Vectors queryHits
#' @importFrom S4Vectors subjectHits
#' @import GenomicRanges
processHits <- function(queryRegions, tx, type) {
  overlaps <- GenomicRanges::findOverlaps(queryRegions, tx)
  overlapsQuery <- queryRegions[S4Vectors::queryHits(overlaps)]
  overlapsTX <- tx[S4Vectors::subjectHits(overlaps)]
  overlapsTX$overlappingQuery <- paste(GenomicRanges::seqnames(overlapsQuery),
                                       GenomicRanges::start(overlapsQuery),
                                       GenomicRanges::end(overlapsQuery),
                                       GenomicRanges::strand(overlapsQuery),
                                       sep=':')
  dt <- data.table::data.table('tx_name' = overlapsTX$tx_name,
                               'query' = overlapsTX$overlappingQuery)
  query <- 'query'
  summary <- dt[,length(unique(query)), by='tx_name']
  colnames(summary) <- c('tx_name', type)
  return(summary)
}


#' getTargetedGenesTable
#'
#' This function provides a list of genes which are targeted by query regions
#' and their corresponding numbers from an input BED file. Then, the hits are
#' categorized by the gene features such as promoters, introns, exons,
#' 5'/3' UTRs and whole transcripts.
#'
#' @param queryRegions GRanges object containing coordinates of input query
#'   regions imported by the \code{\link{importBed}} function
#' @param txdbFeatures A list of GRanges objects where each GRanges object
#'   corresponds to the genomic coordinates of gene features such as promoters,
#'   introns, exons, 5'/3' UTRs and whole transcripts.
#'   This list of GRanges objects are obtained by the function
#'   \code{\link{getTxdbFeaturesFromGRanges}} or \code{\link{getTxdbFeatures}}.
#'
#' @examples
#' data(gff)
#' data(queryRegions)
#' txdbFeatures <- getTxdbFeaturesFromGRanges(gffData = gff)
#' featuresTable <- getTargetedGenesTable(queryRegions = queryRegions,
#'                                        txdbFeatures = txdbFeatures)
#'
#' #or
#' \dontrun{
#' txdb <- GenomicFeatures::makeTxDbFromGRanges(gff)
#' txdbFeatures <- getTxdbFeatures(txdb)
#' featuresTable <- getTargetedGenesTable(queryRegions = queryRegions,
#'                                        txdbFeatures = txdbFeatures)
#'                                        }
#' @return A data.frame object where rows correspond to genes and columns
#'   correspond to gene features
#'
#' @importFrom data.table setkey
#' @importFrom S4Vectors Reduce
#' @export
getTargetedGenesTable <- function (queryRegions, txdbFeatures) {

  tbls <- lapply(X = seq_along(txdbFeatures),
                 FUN = function(i) {
                   processHits(queryRegions = queryRegions,
                               tx = txdbFeatures[[names(txdbFeatures)[i]]],
                               type = names(txdbFeatures)[i])})

  tbls <- lapply(tbls, function(i) data.table::setkey(i, 'tx_name'))
  merged <- S4Vectors::Reduce(function(...) merge(..., all = TRUE), tbls)
  merged[is.na(merged)] <- 0
  return(merged)
}

#' summarizeQueryRegions
#'
#' This function counts number of query regions that overlap with different
#' types of gene features.
#'
#' @param queryRegions GRanges object imported from a BED file using
#'   \code{importBed} function
#' @param txdbFeatures List of GRanges objects - outputs of
#'   \code{getTxdbFeaturesFromGRanges} and \code{getTxdbFeatures} functions
#' @return A data frame with two columns where first column holds features and
#'   second column holds corresponding counts
#' @examples
#' data(gff)
#' data(queryRegions)
#' txdbFeatures <- getTxdbFeaturesFromGRanges(gffData = gff)
#' summary <- summarizeQueryRegions(queryRegions = queryRegions,
#'                                  txdbFeatures = txdbFeatures)
#'
#' @import GenomicRanges
#' @importFrom S4Vectors queryHits
#' @export
summarizeQueryRegions <- function(queryRegions, txdbFeatures) {
  summarize <- function (x) {
    myOverlaps <- GenomicRanges::findOverlaps(queryRegions, x)
    unique(S4Vectors::queryHits(myOverlaps))
  }
  results <- lapply(X = txdbFeatures, FUN = summarize)
  results$NoFeatures <- setdiff(c(1:length(queryRegions)), 
                                as.vector(do.call(c, results)))
  
  df <- t(data.frame(lapply(results, length)))
  colnames(df) <- c('count')
  return(df)
}

#' queryGff
#'
#' This function checks overlaps between the regions in input query and in
#' reference. Input query should be in BED format and reference should be in GFF
#' format. Both data are imported as GRanges object.
#'
#' @param queryRegions GRanges object imported from a BED file using
#'   \code{importBed} function
#' @param gffData GRanges object imported from a GTF file using \code{importGtf}
#'   function
#' @return a GRanges object (a subset of input gff) with an additional column
#'   'overlappingQuery' that contains the coordinates of query regions that
#'   overlap the target annotation features
#'
#' @examples
#' data(queryRegions)
#' data(gff)
#' overlaps <- queryGff(queryRegions = queryRegions, gffData = gff)
#'
#' @import GenomicRanges
#' @importFrom S4Vectors queryHits
#' @importFrom S4Vectors subjectHits
#' @export
queryGff <- function(queryRegions, gffData) {
  #find all overlapping pairs of intervals between gff features and BED file
  overlaps = GenomicRanges::findOverlaps(queryRegions, gffData)
  overlapsQuery = queryRegions[S4Vectors::queryHits(overlaps)]
  overlapsGff = gffData[S4Vectors::subjectHits(overlaps)]
  
  # prepare some data for the overlapping query region 
  # to keep along with overlapping gff features 
  query_mcols <- mcols(overlapsQuery)
  colnames(query_mcols) <- paste0('query_', colnames(query_mcols))
  
  queryData <- cbind(
    data.frame('queryIndex' = S4Vectors::queryHits(overlaps), 
               'queryRange' = paste0(GenomicRanges::seqnames(overlapsQuery), ':', 
                                  GenomicRanges::start(overlapsQuery), '-', 
                                  GenomicRanges::end(overlapsQuery), ':',
                                  GenomicRanges::strand(overlapsQuery))),
    query_mcols)
  
  # update mcols portion of overlapping gff feature data           
  mcols(overlapsGff) <- cbind(mcols(overlapsGff), queryData)
  return (overlapsGff)
}


#' getFeatureBoundaryCoverage
#' 
#' This function extracts the flanking regions of 5' and 3' boundaries of a 
#' given set of genomic features and computes the per-base coverage of query 
#' regions across these boundaries.
#' 
#' @param queryRegions GRanges object imported from a BED file using 
#'   \code{importBed} function
#' @param featureCoords GRanges object containing the target feature coordinates
#' @param flankSize Positive integer that determines the number of base pairs to
#'   extract around a given genomic feature boundary
#' @param boundaryType (Options: fiveprime or threeprime). Denotes which side of
#'   the feature's boundary is to be profiled.
#' @param sampleN A positive integer value less than the total number of featuer
#'   coordinates that determines whether the target feature coordinates should 
#'   be randomly downsampled. If set to 0, no downsampling will happen. If
#' @return a data frame containin three columns. 1. fivePrime: Coverage at 5' 
#'   end of features 2. threePrime: Coverage at 3' end of features; 3. bases: 
#'   distance (in bp) to the boundary
#'
#' @examples
#' data(queryRegions)
#' data(gff)
#' txdb <- GenomicFeatures::makeTxDbFromGRanges(gff)
#' transcriptCoords <- GenomicFeatures::transcripts(txdb)
#' transcriptEndCoverage <- getFeatureBoundaryCoverage (
#'                                      queryRegions = queryRegions,
#'                                     featureCoords = transcriptCoords,
#'                                     flankSize = 100,
#'                                     boundaryType = 'threeprime',
#'                                     sampleN = 1000)
#' @import GenomicRanges
#' @importFrom genomation ScoreMatrix
#' @export
getFeatureBoundaryCoverage <- function (queryRegions,
                                         featureCoords,
                                         flankSize = 500,
                                         boundaryType, 
                                         sampleN = 0) {
  
  if (sampleN > 0 && sampleN < length(featureCoords)) {
    featureCoords <- sort(featureCoords[sample(length(featureCoords), sampleN)])
  }
  
  if (boundaryType == 'fiveprime') {
    flanks <- GenomicRanges::flank(x = featureCoords,
                                   width = flankSize,
                                   start = TRUE,
                                   both = TRUE)
  } else if (boundaryType == 'threeprime') { 
    flanks <- GenomicRanges::flank(x = featureCoords,
                                   width = flankSize,
                                   start = FALSE,
                                   both = TRUE)
  } else {
    stop ("please indicate either threeprime or fiveprime for boundary type\n")
  }
  
  sm <- genomation::ScoreMatrix(target = queryRegions,
                                windows = flanks,
                                strand.aware = TRUE)
  
  return (data.frame('bases' = c(-flankSize:(flankSize-1)), 
                     'meanCoverage' = colMeans(sm), 
                     'standardError' = apply(sm, 2, plotrix::std.error))
  )
}
  

#' getFeatureBoundaryCoverageBin
#'
#' This function extracts the flanking regions of 5' and 3' boundaries of a
#' given set of genomic features, splits them into 100 equally sized bins and
#' computes the per-bin coverage of query regions across these boundaries.
#'
#' @param queryRegions GRanges object imported from a BED file using
#'   \code{importBed} function
#' @param featureCoords GRanges object containing the target feature coordinates
#' @param flankSize Positive integer that determines the number of base pairs to
#'   extract around a given genomic feature boundary
#' @param sampleN A positive integer value less than the total number of featuer
#'   coordinates that determines whether the target feature coordinates should
#'   be randomly downsampled. If set to 0, no downsampling will happen. If
#' @return a data frame containin three columns. 1. fivePrime: Coverage at 5'
#'   end of features 2. threePrime: Coverage at 3' end of features; 3. bases:
#'   distance (in bp) to the boundary
#'
#' @examples
#' data(queryRegions)
#' data(gff)
#' txdb <- GenomicFeatures::makeTxDbFromGRanges(gff)
#' transcriptCoords <- GenomicFeatures::transcripts(txdb)
#' transcriptEndCoverageBin <- getFeatureBoundaryCoverageBin (
#'                                      queryRegions = queryRegions,
#'                                     featureCoords = transcriptCoords,
#'                                     flankSize = 100,
#'                                     sampleN = 1000)
#' @import GenomicRanges
#' @importFrom genomation ScoreMatrix
#' @export
getFeatureBoundaryCoverageBin <- function (queryRegions,
                                           featureCoords,
                                           flankSize = 50,
                                           sampleN = 0) {

  if (sampleN > 0 && sampleN < length(featureCoords)) {
    featureCoords <- sort(featureCoords[sample(length(featureCoords), sampleN)])
  }

  #flanking regions at/around 5' site of the features
  fivePrimeFlanks <- GenomicRanges::flank(x = featureCoords,
                                          width = flankSize,
                                          start = TRUE,
                                          both = TRUE
  )
  #flanking regions at/around 3' site of the features
  threePrimeFlanks <- GenomicRanges::flank(x = featureCoords,
                                           width = flankSize,
                                           start = FALSE,
                                           both = TRUE)

  cvgFivePrime <- genomation::ScoreMatrixBin(target = queryRegions,
                                             windows = fivePrimeFlanks,
                                             bin.num = 100,
                                             bin.op = 'max',
                                             strand.aware = TRUE)

  cvgThreePrime <- genomation::ScoreMatrixBin(target = queryRegions,
                                              windows = threePrimeFlanks,
                                              bin.num = 100,
                                              bin.op = 'max',
                                              strand.aware = TRUE)

  mdata <- data.frame('fivePrime' = colSums(cvgFivePrime),
                      'threePrime' = colSums(cvgThreePrime),
                      'bins' =  c(-50:-1, 1:50))

  return(mdata)
}


#' calculateCoverageProfile
#' 
#' This function checks overlaps between input query regions and annotation 
#' features, and then calculates coverage profile along target regions.
#' 
#' @param queryRegions GRanges object imported from a BED file using 
#'   \code{importBed} function
#' @param targetRegions GRanges object containing genomic coordinates of a 
#'   target feature (e.g. exons)
#' @param sampleN If set to a positive integer, \code{targetRegions} will be 
#'   downsampled to \code{sampleN} regions
#' @param bin.num Positive integer value (default: 100) to determine how many
#'   bins the targetRegions should be split into (See
#'   genomation::ScoreMatrixBin)
#' @param bin.op The operation to apply for each bin: 'min', 'max', or 'mean'
#'   (default: mean). (See genomation::ScoreMatrixBin)
#' @param strand.aware TRUE/FALSE (default: TRUE) The strands of target regions
#'   are considered.
#' @return A ScoreMatrix object returned by \code{genomation::ScoreMatrixBin} 
#'   function. Target regions are divided into 100 equal sized bins and coverage
#'   level is calculated in a strand-specific manner.
#' @examples
#' data(gff)
#' data(queryRegions)
#' txdbFeatures <- getTxdbFeaturesFromGRanges(gffData = gff)
#' df <- calculateCoverageProfile(queryRegions = queryRegions,
#'                               targetRegions = txdbFeatures$exons,
#'                                     sampleN = 1000)
#' @importFrom genomation ScoreMatrixBin
#' @import GenomicRanges
#' @export
calculateCoverageProfile = function (queryRegions, 
                                     targetRegions, 
                                     sampleN = 0, 
                                     bin.num = 100, 
                                     bin.op = 'mean', 
                                     strand.aware = TRUE){
  #remove windows shorter than bin.num
  #because the windows have to be separated into at least 'bin.num' bins
  windows <- targetRegions[GenomicRanges::width(targetRegions) >= bin.num]

  if (length(windows) > 0) {
    if (sampleN > 0 && sampleN < length(windows)) {
      windows <- sort(windows[sample(length(windows), sampleN)])
    }
    sm <- genomation::ScoreMatrixBin(target = queryRegions,
                                     windows = windows,
                                     bin.num = bin.num,
                                     bin.op = bin.op,
                                     strand.aware = strand.aware)
    return(sm)
    } else {
    stop("Cannot compute coverage profile for target regions.\n
         There are no target regions longer than",bin.num,"\n")
  }
}

#' calculateCoverageProfileList
#' 
#' This function checks overlaps between input query regions and a target list 
#' of annotation features, and then calculates the coverage profile along the 
#' target regions.
#' 
#' @param queryRegions GRanges object imported from a BED file using 
#'   \code{importBed} function
#' @param targetRegionsList A list of GRanges objects containing genomic 
#'   coordinates of target features (e.g. transcripts, exons, introns)
#' @param sampleN If set to a positive integer, \code{targetRegions} will be 
#'   downsampled to \code{sampleN} regions
#' @param bin.num Positive integer value (default: 100) to determine how many
#'   bins the targetRegions should be split into (See
#'   genomation::ScoreMatrixBin)
#' @param bin.op The operation to apply for each bin: 'min', 'max', or 'mean'
#'   (default: mean). (See genomation::ScoreMatrixBin)
#' @param strand.aware TRUE/FALSE (default: TRUE) The strands of target regions
#'   are considered.
#' @return A data.frame consisting of four columns: 1. bins level 2.
#'   meanCoverage 3. standardError 4. feature Target regions are divided into
#'   100 equal sized bins and coverage level is summarized in a strand-specific
#'   manner using the \code{genomation::ScoreMatrixBin} function. For each bin,
#'   mean coverage score and the standard error of the mean coverage score is
#'   calculated (\code{plotrix::std.error})
#' @importFrom plotrix std.error
#' @examples
#' data(gff)
#' data(queryRegions)
#' txdbFeatures <- getTxdbFeaturesFromGRanges(gffData = gff)
#' dfList <- calculateCoverageProfileList(queryRegions = queryRegions,
#'                               targetRegionsList = txdbFeatures,
#'                                     sampleN = 1000)
#' @export
calculateCoverageProfileList <- function (queryRegions,
                                          targetRegionsList,
                                          sampleN = 0, 
                                          bin.num = 100,
                                          bin.op = 'mean',
                                          strand.aware = TRUE) {
  results <- lapply(X = names(targetRegionsList),
                    FUN = function(x) { 
                      sm <- calculateCoverageProfile(queryRegions = queryRegions,
                                                     targetRegions = targetRegionsList[[x]],
                                                     sampleN = sampleN, 
                                                     bin.num = bin.num, 
                                                     bin.op = bin.op, 
                                                     strand.aware = strand.aware)
                      mdata <- data.frame('bins' = c(1:bin.num), 
                                          'meanCoverage' = colMeans(sm), 
                                          'standardError' = apply(sm, 2, plotrix::std.error),
                                          'feature' = x)
                    })
  return(do.call(rbind, results))
}

#' calculateCoverageProfileListFromTxdb
#'
#' This function is deprecated. Use ?calculateCoverageProfileList instead. 
#' 
#' @export
calculateCoverageProfileListFromTxdb <- function () {
  .Deprecated('calculateCoverageProfileList')
}

#' calculateCoverageProfileFromTxdb
#'
#' This function is deprecated. Use ?calculateCoverageProfile instead. 
#' 
#' @export
calculateCoverageProfileFromTxdb <- function () {
  .Deprecated('calculateCoverageProfile')
}

#' checkSeqDb
#' 
#' Given a string that denotes a genome version (e.g. hg19)
#' returns the BSgenome object matching the genome version 
#' that are available in BSgenome::available.genomes() 
#' 
#' @param genomeVersion String that denotes genome version. 
#' To unambigously select a BSgenome object, provide a string
#' that matches the end of the available genomes at: 
#' BSgenome::available.genomes().
#' @return Returns a BSgenome object that uniquely matches the 
#' genomeVersion. 
#' @examples 
#' checkSeqDb('hg19')
#' @importFrom BSgenome available.genomes
#' @export
checkSeqDb <- function(genomeVersion) {
  availableGenomes <- BSgenome::available.genomes()
  db <- grep(paste0(genomeVersion, '$'),
         BSgenome::available.genomes(),
         value = T)
  if (length(db) == 0) {
    stop(
      "Can't find a genome sequence with the given genome version:",
      genomeVersion,
      "\n\tRun BSgenome::available.genomes() to see which ones are available"
    )
  } else if (length(db) > 1) {
    stop(
      "Can't match the genome version to an unambigous genome sequence object.",
      "\nmatching genomes are: ",
      paste(db, collapse = '\t'),
      "\nPlease update your genomeVersion input to uniquely match",
      "\nAlso see BSgenome::available.genomes()"
    )
  }
  
  if (!requireNamespace(db, quietly = TRUE)) {
    stop(
      "Can't find the package ",
      db,
      " in your system. \n Please install via => ",
      paste0("BiocManager::install('", db, "')")
    )
  }
  
  message(date(), " => Returning ", db, " as BSgenome object")
  
  require(db, character.only = TRUE)
  return(get(db))
}

#' Generate a RCAS Report for a list of transcriptome-level segments
#' 
#' This is the main report generation function for RCAS. This function takes a
#' BED file, a GTF file to create a summary report regarding the annotation data
#' that overlap the input BED file, enrichment analysis for functional terms,
#' and motif analysis.
#' 
#' @param queryFilePath a BED format file which contains genomic coordinates of 
#'   protein-RNA binding sites
#' @param gffFilePath A GTF format file which contains genome annotations 
#'   (preferably from ENSEMBL)
#' @param annotationSummary TRUE/FALSE (default: TRUE) A switch to decide if 
#'   RCAS should provide annotation summaries from overlap operations
#' @param goAnalysis TRUE/FALSE (default: TRUE) A switch to decide if RCAS 
#'   should run GO term enrichment analysis
#' @param motifAnalysis TRUE/FALSE (default: TRUE) A switch to decide if RCAS 
#'   should run motif analysis
#' @param genomeVersion  A character string to denote for which genome version 
#'   the analysis is being done. 
#' @param outDir Path to the output directory. (default: current working 
#'   directory)
#' @param printProcessedTables boolean value (default: FALSE). If set to TRUE, 
#'   raw data tables that are used for plots/tables will be printed to text
#'   files.
#' @param sampleN integer value (default: 0). A parameter to determine if the 
#'   input query regions should be downsampled to a smaller size in order to 
#'   make report generation quicker. When set to 0, downsampling won't be done. 
#'   To activate the sampling a positive integer value that is smaller than the
#'   total number of query regions should be given.
#' @param quiet boolean value (default: FALSE). If set to TRUE, progress bars 
#' and chunk labels will be suppressed while knitting the Rmd file. 
#' @param selfContained boolean value (default: TRUE). By default, the generated
#' html file will be self-contained, which means that all figures and tables 
#' will be embedded in a single html file with no external dependencies 
#' (See rmarkdown::html_document) 
#' @return An html generated using rmarkdown/knitr/pandoc that contains 
#'   interactive figures, tables, and text that provide an overview of the 
#'   experiment
#' @examples
#' #Default run will generate a report using built-in test data for hg19 genome.
#' \dontrun{
#' runReport()
#' }
#' 
#' #A custom run for human
#' \dontrun{
#' runReport( queryFilePath = 'input.BED',
#'            gffFilePath = 'annotation.gtf', 
#'            genomeVersion = 'hg19')
#'            }
#' # To turn off certain modules of the report
#' \dontrun{
#' runReport( queryFilePath = 'input.BED',
#'            gffFilePath = 'annotation.gtf',
#'            motifAnalysis = FALSE,
#'            goAnalysis = FALSE )
#'            }
#' @import rmarkdown
#' @import knitr
#' @importFrom plotly plot_ly
#' @importFrom plotly add_trace
#' @importFrom plotly add_ribbons
#' @importFrom plotly layout
#' @import DT
#' @export
runReport <- function(queryFilePath = 'testdata',
                      gffFilePath = 'testdata',
                      annotationSummary = TRUE,
                      goAnalysis = TRUE,
                      motifAnalysis = TRUE,
                      genomeVersion = 'hg19',
                      outDir = getwd(),
                      printProcessedTables = FALSE,
                      sampleN = 0,
                      quiet = FALSE,
                      selfContained = TRUE) {

  db <- checkSeqDb(genomeVersion)
  # get species name 
  # this is needed for gprofiler functional enrichment 
  fields <- unlist(strsplit(db@metadata$organism, ' '))
  species <- tolower(paste0(unlist(strsplit(fields[1], ''))[1], 
                 fields[2]))
  
  if(queryFilePath != 'testdata') {
    queryFilePath <- normalizePath(queryFilePath)
  }

  if(gffFilePath != 'testdata') {
    gffFilePath <- normalizePath(gffFilePath)
  }
  
  reportFile <- system.file("reporting_scripts", "report.Rmd", package='RCAS')
  headerFile <- system.file("reporting_scripts", "header.html", package='RCAS')
  footerFile <- system.file("reporting_scripts", "footer.html", package='RCAS')
  
  outFile <- paste0(basename(queryFilePath), '.RCAS.report.html')

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
    params = list(query = queryFilePath,
                  gff = gffFilePath,
                  annotationSummary = annotationSummary,
                  goAnalysis = goAnalysis,
                  motifAnalysis = motifAnalysis,
                  genomeVersion = genomeVersion,
                  species = species,
                  printProcessedTables = printProcessedTables,
                  sampleN = sampleN,
                  workdir = outDir),
    quiet = quiet
    )
}

