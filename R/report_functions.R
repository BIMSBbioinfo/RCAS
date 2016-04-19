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
#'
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
                       overwriteObjectAsRds = FALSE, keepStandardChr = TRUE) {

  rdsFilePath <- paste0(filePath, ".granges.rds")

  if (readFromRds == TRUE && file.exists(rdsFilePath)){
    cat('Reading existing granges.rds object from',rdsFilePath,'\n')
    gff <- readRDS(rdsFilePath)
  } else {
    cat('importing gtf file from', filePath,'\n')
    gff <- rtracklayer::import.gff(filePath)
  }

  if (!is.null(gff)) {
    if (keepStandardChr == TRUE) {
      cat('Keeping standard chromosomes only\n')
      GenomeInfoDb::seqlevelsStyle(gff) = 'UCSC'
      gff = GenomeInfoDb::keepStandardChromosomes(gff)
    }

    if (saveObjectAsRds == TRUE) {
      if (!file.exists(rdsFilePath)){
        cat('Saving gff object in RDS file at',rdsFilePath,'\n')
        saveRDS(gff, file=rdsFilePath)
      } else {
        if (overwriteObjectAsRds == TRUE) {
          cat('Overwriting gff object in RDS file at',rdsFilePath,'\n')
          saveRDS(gff, file=rdsFilePath)
        } else {
          message('File ',rdsFilePath,' already exists.\n
                  Use overwriteObjectAsRds = TRUE to overwrite the file\n')
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
#' @param filePath Path to a GTF file
#' @param sampleN A positive integer value. The number of intervals in the input
#'   BED file are randomly downsampled to include intervals as many as
#'   \code{sampleN}. The input will be downsampled only if this value is larger
#'   than zero and less than the total number of input intervals.
#' @param keepStandardChr TRUE/FALSE (default:TRUE). If set to TRUE, will
#'   convert the \code{seqlevelsStyle} to 'UCSC' and apply
#'   \code{keepStandardChromosomes} function to only keep data from the standard
#'   chromosomes
#'
#' @return A \code{GRanges} object containing the coordinates of the intervals
#'   from an input BED file
#'
#' @examples
#' input <- system.file('testfile.bed', package='RCAS')
#' importBed(filePath = input, keepStandardChr = TRUE)
#'
#' @importFrom rtracklayer import.bed
#' @importFrom GenomeInfoDb seqlevelsStyle
#' @importFrom GenomeInfoDb keepStandardChromosomes
#' @export
importBed <- function (filePath, sampleN = 0, keepStandardChr = TRUE) {

  if (file.exists(filePath)) {
    data = rtracklayer::import.bed(filePath)
    if (keepStandardChr == TRUE) {
      GenomeInfoDb::seqlevelsStyle(data) <- 'UCSC'
      cat('Keeping standard chromosomes only\n')
      data <- GenomeInfoDb::keepStandardChromosomes(data)
    }

    if (sampleN > 0 && sampleN < length(data)) {
      cat ('Downsampling intervals to size:',sampleN,'\n')
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
#' This function takes as input a txdb object from GenomicFeatures library. Then
#' extracts the coordinates of gene features such as promoters, introns, exons,
#' intron-exon boundaries, 5'/3' UTRs, and whole transcripts.
#'
#' @param txdb A txdb object imported by GenomicFeatures::makeTxDb family of
#'   functions
#'
#' @examples
#' data(gff)
#' txdb <- GenomicFeatures::makeTxDbFromGRanges(gff)
#' txdbFeatures <- getTxdbFeatures(txdb)
#'
#' @return A list of GRanges objects
#'
#' @import GenomicFeatures
#' @import GenomicRanges
#' @importFrom BiocGenerics unlist
#' @export
getTxdbFeatures <- function (txdb) {

  transcripts <- GenomicFeatures::transcripts(txdb)

  tmp <- GenomicFeatures::exonsBy(x = txdb, by = "tx", use.names = TRUE)
  exons <- BiocGenerics::unlist(tmp)
  exons$tx_name <- names(exons)

  tmp <- GenomicFeatures::intronsByTranscript(x = txdb, use.names = TRUE)
  introns <- BiocGenerics::unlist(tmp)
  introns$tx_name <- names(introns)

  tmp1 <- GenomicRanges::flank(x = introns[GenomicRanges::width(introns) > 99],
                               width = 50,
                               start = TRUE,
                               both = TRUE)
  tmp2 <- GenomicRanges::flank(x = introns[GenomicRanges::width(introns) > 99],
                               width=50,
                               start=FALSE,
                               both=TRUE)
  exonIntronBoundaries <- c(tmp1, tmp2)
  exonIntronBoundaries$tx_name <- names(exonIntronBoundaries)

  promoters <- GenomicFeatures::promoters(txdb)

  tmp <- GenomicFeatures::fiveUTRsByTranscript(x = txdb, use.names = TRUE)
  fiveUTRs <- BiocGenerics::unlist(tmp)
  fiveUTRs$tx_name <- names(fiveUTRs)

  tmp <- GenomicFeatures::threeUTRsByTranscript(x = txdb, use.names = TRUE)
  threeUTRs <- BiocGenerics::unlist(tmp)
  threeUTRs$tx_name <- names(threeUTRs)

  tmp <- GenomicFeatures::cdsBy(x = txdb, by = "tx", use.names = TRUE)
  cds <- BiocGenerics::unlist(tmp)
  cds$tx_name <- names(cds)

  txdbFeatures <- list(
    'transcripts' = transcripts,
    'exons'       = exons,
    'promoters'   = promoters,
    'fiveUTRs'    = fiveUTRs,
    'introns'     = introns,
    'exonIntronBoundaries' = exonIntronBoundaries,
    'cds'         = cds,
    'threeUTRs'   = threeUTRs
  )
  return(txdbFeatures)
}

#' getTxdbFeaturesFromGff
#'
#' This function takes as input a GRanges object from the output of
#' \code{importGtf} function. Then extracts the coordinates of gene features
#' such as promoters, introns, exons, intron-exon boundaries, 5'/3' UTRs and
#' whole transcripts.
#'
#' @param gff A GRanges object imported by \code{importGtf} function
#'
#' @examples
#' data(gff)
#' txdbFeatures <- getTxdbFeaturesFromGff(gff)
#'
#' @return A list of GRanges objects
#'
#' @import GenomicFeatures
#' @import GenomicRanges
#' @importFrom BiocGenerics unlist
#' @export
getTxdbFeaturesFromGff <- function (gff) {

  txdb <- GenomicFeatures::makeTxDbFromGRanges(gff)

  transcripts <- GenomicFeatures::transcripts(txdb)
  m <- match(transcripts$tx_name, gff$transcript_id)
  transcripts$gene_name <- gff[m]$gene_name

  tmp <- GenomicFeatures::exonsBy(x = txdb, by = "tx", use.names = TRUE)
  exons <- BiocGenerics::unlist(tmp)
  exons$tx_name <- names(exons)
  exons$gene_name <- gff[match(names(exons), gff$transcript_id)]$gene_name

  tmp <- GenomicFeatures::intronsByTranscript(x = txdb, use.names = TRUE)
  introns <- BiocGenerics::unlist(tmp)
  introns$tx_name <- names(introns)
  m <- match(names(introns), gff$transcript_id)
  introns$gene_name <- gff[m]$gene_name

  tmp1 <- GenomicRanges::flank(introns[GenomicRanges::width(introns) >= 100],
                               width = 50, start = TRUE, both = TRUE)
  tmp2 <- GenomicRanges::flank(introns[GenomicRanges::width(introns) >= 100],
                               width = 50, start = FALSE, both = TRUE)
  exonIntronBoundaries <- c(tmp1, tmp2)

  promoters <- GenomicFeatures::promoters(txdb)
  m <- match(promoters$tx_name, gff$transcript_id)
  promoters$gene_name <- gff[m]$gene_name

  tmp <- GenomicFeatures::fiveUTRsByTranscript(x = txdb, use.names = TRUE)
  fiveUTRs <- BiocGenerics::unlist(tmp)
  fiveUTRs$tx_name <- names(fiveUTRs)
  m <- match(names(fiveUTRs), gff$transcript_id)
  fiveUTRs$gene_name <- gff[m]$gene_name

  tmp <- GenomicFeatures::threeUTRsByTranscript(x = txdb, use.names = TRUE)
  threeUTRs <- BiocGenerics::unlist(tmp)
  threeUTRs$tx_name <- names(threeUTRs)
  m <- match(names(threeUTRs), gff$transcript_id)
  threeUTRs$gene_name <- gff[m]$gene_name

  tmp <- GenomicFeatures::cdsBy(x = txdb, by = "tx", use.names = TRUE)
  cds <- BiocGenerics::unlist(tmp)
  cds$tx_name <- names(cds)
  cds$gene_name <- gff[match(names(cds), gff$transcript_id)]$gene_name

  txdbFeatures = list(
    'transcripts' = transcripts,
    'exons'       = exons,
    'promoters'   = promoters,
    'fiveUTRs'    = fiveUTRs,
    'introns'     = introns,
    'exonIntronBoundaries' = exonIntronBoundaries,
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
  summary <- dt[,length(unique(query)), by='tx_name']
  colnames(summary) <- c('tx_name', type)
  return(summary)
}


#' getTargetedGenesTable
#'
#' This function provides a list of genes which are targeted by query regions
#' and their corresponding numbers from an input BED file. Then, the hits are
#' categorized by the gene features such as promoters, introns, exons,
#' intron-exon boundaries, 5'/3' UTRs and whole transcripts.
#'
#' @param queryRegions GRanges object containing coordinates of input query
#'   regions imported by the \code{\link{importBed}} function
#' @param txdbFeatures A list of GRanges objects where each GRanges object
#'   corresponds to the genomic coordinates of gene features such as promoters,
#'   introns, exons, intron-exon boundaries, 5'/3' UTRs and whole transcripts.
#'   This list of GRanges objects are obtained by the function
#'   \code{\link{getTxdbFeaturesFromGff}} or \code{\link{getTxdbFeatures}}.
#'
#' @examples
#' data(gff)
#' data(queryRegions)
#' txdbFeatures <- getTxdbFeaturesFromGff(gff)
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
#'   \code{getTxdbFeaturesFromGff} and \code{getTxdbFeatures} functions
#' @return A data frame with two columns where first column holds features and
#'   second column holds corresponding counts
#' @examples
#' data(gff)
#' data(queryRegions)
#' txdbFeatures <- getTxdbFeaturesFromGff(gff)
#' summary <- summarizeQueryRegions(queryRegions = queryRegions,
#'                                  txdbFeatures = txdbFeatures)
#'
#' @import GenomicRanges
#' @importFrom S4Vectors queryHits
#' @export
summarizeQueryRegions <- function(queryRegions, txdbFeatures) {
  summarize <- function (x) {
    myOverlaps <- GenomicRanges::findOverlaps(queryRegions, x)
    tmp <- S4Vectors::queryHits(myOverlaps)
    return(length(unique(tmp)))
  }
  results = lapply(X = txdbFeatures, FUN = summarize)
  df <- t(data.frame(results))
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
#' @param gff GRanges object imported from a GTF file using \code{importGtf}
#'   function
#' @return a GRanges object (a subset of input gff) with an additional column
#'   'overlappingQuery' that contains the coordinates of query regions that
#'   overlap the target annotation features
#'
#' @examples
#' data(queryRegions)
#' data(gff)
#' overlaps <- queryGff(queryRegions = queryRegions, gff = gff)
#'
#' @import GenomicRanges
#' @importFrom S4Vectors queryHits
#' @importFrom S4Vectors subjectHits
#' @export
queryGff <- function(queryRegions, gff) {
  #find all overlapping pairs of intervals between gff features and BED file
  overlaps = GenomicRanges::findOverlaps(queryRegions, gff)
  overlapsQuery = queryRegions[S4Vectors::queryHits(overlaps)]
  overlapsGff = gff[S4Vectors::subjectHits(overlaps)]
  overlapsGff$overlappingQuery = paste(GenomicRanges::seqnames(overlapsQuery),
                                       GenomicRanges::start(overlapsQuery),
                                       GenomicRanges::end(overlapsQuery),
                                       GenomicRanges::strand(overlapsQuery),
                                       sep=':')
  return (overlapsGff)
}

#' calculateCoverageProfile
#'
#' This function checks overlaps betwenn input query regions and annotation
#' features, and then calculates coverage profile along target regions.
#'
#' @param queryRegions GRanges object imported from a BED file using
#'   \code{importBed} function
#' @param targetRegions GRanges object containing genomic coordinates of a
#'   target feature (e.g. exons)
#' @param sampleN If set to a positive integer, \code{targetRegions} will be
#'   downsampled to \code{sampleN} regions
#'
#' @return A data.frame object consisting of two columns: 1. coverage level 2.
#'   bins. Target regions are divided into 100 equal sized bins and coverage
#'   level is summarized in a strand-specific manner using the
#'   \code{genomation::ScoreMatrixBin} function.
#' @examples
#' data(gff)
#' data(queryRegions)
#' txdbFeatures <- getTxdbFeaturesFromGff(gff)
#' df <- calculateCoverageProfile(queryRegions = queryRegions,
#'                               targetRegions = txdbFeatures$exons,
#'                                     sampleN = 1000)
#'
#' @importFrom genomation ScoreMatrixBin
#' @import GenomicRanges
#' @export
calculateCoverageProfile = function (queryRegions, targetRegions, sampleN = 0){
  #remove windows shorter than 100 bp
  #because the windows have to be separated into at least 100 bins
  windows <- targetRegions[GenomicRanges::width(targetRegions) >= 100]

  if (length(windows) > 0) {
    if (sampleN > 0 && sampleN < length(windows)) {
      windows <- windows[sample(length(windows), sampleN)]
    }
    sm <- genomation::ScoreMatrixBin(target = queryRegions,
                                     windows = windows,
                                     bin.num = 100,
                                     strand.aware = TRUE)
    mdata <- as.data.frame(colSums(sm))
    mdata$bins <- c(1:100)
    colnames(mdata) <- c('coverage', 'bins')
    return(mdata)
    } else {
    stop("Cannot compute coverage profile for target regions.\n
         There are no target regions longer than 100 bp\n")
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
#'
#' @return A list of data.frame objects consisting of two columns: 1. coverage
#'   level 2. bins. Target regions are divided into 100 equal sized bins and
#'   coverage level is summarized in a strand-specific manner using the
#'   \code{genomation::ScoreMatrixBin} function.
#' @examples
#' data(gff)
#' data(queryRegions)
#' txdbFeatures <- getTxdbFeaturesFromGff(gff)
#' list.df <- calculateCoverageProfileList(queryRegions = queryRegions,
#'                               targetRegionsList = txdbFeatures,
#'                                     sampleN = 1000)
#' @export
calculateCoverageProfileList <- function (queryRegions,
                                          targetRegionsList,
                                          sampleN = 0) {
  lapply(X = targetRegionsList,
       FUN = function(x) { calculateCoverageProfile(queryRegions = queryRegions,
                                                    targetRegions = x,
                                                    sampleN = sampleN) })
}

#' calculateCoverageProfileListFromTxdb
#'
#' This function overlaps the input query regions with a target list of
#' annotation features and calculates the coverage profile along the target
#' regions.
#'
#' @param queryRegions GRanges object imported from a BED file using
#'   \code{importBed} function
#' @param txdb A txdb object obtained by using \code{GenomicFeatures::makeTxDb}
#'   family of functions
#' @param sampleN If set to a positive integer, \code{targetRegions} will be
#'   downsampled to \code{sampleN} regions
#'
#' @return A list of data.frame objects consisting of two columns: 1. coverage
#'   level 2. bins. The target regions are divided into 100 equal sized bins and
#'   coverage level is summarized in a strand-specific manner using the
#'   \code{genomation::ScoreMatrixBin} function.
#' @examples
#' data(gff)
#' data(queryRegions)
#' txdb <- GenomicFeatures::makeTxDbFromGRanges(gff)
#' df <- calculateCoverageProfileListFromTxdb(queryRegions = queryRegions,
#'                                                    txdb = txdb,
#'                                                 sampleN = 1000)
#'
#' @export
calculateCoverageProfileListFromTxdb <- function (queryRegions,
                                                  txdb,
                                                  sampleN = 0) {
  txdbFeatures = getTxdbFeatures(txdb = txdb)
  lapply(X = txdbFeatures,
       FUN = function(x) {
         calculateCoverageProfile(queryRegions,
                                  x,
                                  sampleN = sampleN) })
}

#' calculateCoverageProfileFromTxdb
#'
#' This function overlaps the input query regions with a target list of
#' annotation features and calculates the coverage profile along the target
#' regions.
#'
#' @param queryRegions GRanges object imported from a BED file using
#'   \code{importBed} function
#' @param txdb A txdb object obtained by using \code{GenomicFeatures::makeTxDb}
#'   family of functions
#' @param sampleN If set to a positive integer, the targetRegions will be
#'   downsampled to \code{sampleN} regions
#' @param type A character string defining the type of gene feature for which a
#'   profile should be calculated. The options are: transcripts, exons, introns,
#'   exonIntronBoundaries, promoters, fiveUTRs, threeUTRs, and cds.
#' @return A data.frame object consisting of two columns: 1. coverage level 2.
#'   bins. The target regions are divided into 100 equal sized bins and coverage
#'   level is summarized in a strand-specific manner using the
#'   \code{genomation::ScoreMatrixBin} function.
#' @examples
#' data(gff)
#' data(queryRegions)
#' txdb <- GenomicFeatures::makeTxDbFromGRanges(gff)
#' df <- calculateCoverageProfileFromTxdb(queryRegions = queryRegions,
#'                                                type = 'exons',
#'                                                txdb = txdb,
#'                                             sampleN = 1000)
#' @import GenomicFeatures
#' @import GenomicRanges
#' @importFrom BiocGenerics unlist
#' @export
calculateCoverageProfileFromTxdb <- function (queryRegions,
                                              txdb,
                                              type,
                                              sampleN = 0) {

  if (type == 'transcripts') {
    targetRegions <- GenomicFeatures::transcripts(txdb)
  } else if (type == 'exons') {
    tmp <- GenomicFeatures::exonsBy(x = txdb, by = "tx", use.names = TRUE)
    targetRegions <- BiocGenerics::unlist(tmp)
  } else if (type == 'introns') {
    tmp <- GenomicFeatures::intronsByTranscript(x = txdb, use.names = TRUE)
    targetRegions <- BiocGenerics::unlist(tmp)
  } else if (type == 'exonIntronBoundaries') {
    tmp <- GenomicFeatures::intronsByTranscript(x = txdb, use.names = TRUE)
    introns <- BiocGenerics::unlist(tmp)
    tmp1 <- GenomicRanges::flank(introns[GenomicRanges::width(introns) >= 100],
                                 width = 50,
                                 start = TRUE,
                                 both = TRUE)
    tmp2 <- GenomicRanges::flank(introns[GenomicRanges::width(introns) >= 100],
                                 width = 50,
                                 start = FALSE,
                                 both = TRUE)
    targetRegions <- c(tmp1, tmp2)
  } else if (type == 'promoters') {
    targetRegions <- GenomicFeatures::promoters(txdb)
  } else if (type == 'fiveUTRs') {
    tmp <- GenomicFeatures::fiveUTRsByTranscript(x = txdb, use.names = TRUE)
    targetRegions <- BiocGenerics::unlist(tmp)
  } else if (type == 'threeUTRs') {
    tmp <- GenomicFeatures::threeUTRsByTranscript(x = txdb, use.names = TRUE)
    targetRegions <- BiocGenerics::unlist()
  } else if (type == 'cds') {
    tmp <- GenomicFeatures::cdsBy(x = txdb, by = "tx", use.names = TRUE)
    targetRegions <- BiocGenerics::unlist(tmp)
  } else {
    stop ("Can calculate coverage profiles for only:
          transcripts, exons, introns, exonIntronBoundaries,
          promoters, fiveUTRs, threeUTRs, and cds")
  }
  result = calculateCoverageProfile(queryRegions = queryRegions,
                                    targetRegions = targetRegions,
                                    sampleN = sampleN)
  return (result)
}

findLongLines <- function (myfile, lineLimit = 80) {
  counts <- lapply(X = readLines(myfile), FUN = function (x) {l=nchar(x)})
  names(counts) <- c(1:length(counts))
  return(names(counts)[counts > lineLimit])
}

#' Generate a RCAS Report for a CLIP-Seq-based experiment
#'
#' This is the main report generation function for RCAS. This function can take
#' a BED file, a GTF file and optionally an MSIGDB annotation; and use these
#' input to run multiple RCAS functions to create a summary report regarding the
#' annotation data that overlap the input BED file, enrichment analysis for GO
#' terms, gene sets from MSIGDB, and motif analysis.
#'
#' @param queryFilePath a BED format file which contains genomic coordinates of
#'   protein-RNA binding sites
#' @param gffFilePath A GTF format file which contains genome annotations
#'   (preferably from ENSEMBL)
#' @param msigdbFilePath Gene set annotations for Homo sapiens from Molecular
#'   Signatures Database
#' @param annotationSummary TRUE/FALSE (default: TRUE) A switch to decide if
#'   RCAS should provide annotation summaries from overlap operations
#' @param goAnalysis TRUE/FALSE (default: TRUE) A switch to decide if RCAS
#'   should run GO term enrichment analysis
#' @param msigdbAnalysis TRUE/FALSE (default: TRUE) A switch to decide if RCAS
#'   should run gene set enrichment analysis
#' @param motifAnalysis TRUE/FALSE (default: TRUE) A switch to decide if RCAS
#'   should run motif analysis
#' @param genomeVersion  A character string to denote for which genome version
#'   the analysis is being done. Available options are hg19 (human), mm9
#'   (mouse), ce10 (worm) and dm3 (fly).
#' @param outDir Path to the output directory. (default: current working
#'   directory)
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
#'            msigdbFilePath = 'human_msigdb.gmt')
#'            }
#' # To turn off certain modules of the report
#' \dontrun{
#' runReport( queryFilePath = 'input.BED',
#'            gffFilePath = 'annotation.gtf',
#'            msigdbFilePath = 'human_msigdb.gmt',
#'            motifAnalysis = FALSE,
#'            goAnalysis = FALSE )
#'            }
#' # To run the pipeline for species other than human
#' # If the msigdb module is needed, the msigdbFilePath
#' # must be set to the MSIGDB annotations for 'human'.
#' # MSIGDB datasets for other species will be calculated
#' # in the background using the createOrthologousMsigdbDataset
#' # function
#' \dontrun{
#' runReport( queryFilePath = 'input.mm9.BED',
#'            gffFilePath = 'annotation.mm9.gtf',
#'            msigdbFilePath = 'msigdb.human.gmt',
#'            genomeVersion = 'mm9' )
#'            }
#' @import rmarkdown
#' @import plotly
#' @import DT
#' @export
runReport <- function(queryFilePath = 'testdata',
                      gffFilePath = 'testdata',
                      msigdbFilePath = 'testdata',
                      annotationSummary = TRUE,
                      goAnalysis = TRUE,
                      msigdbAnalysis = TRUE,
                      motifAnalysis = TRUE,
                      genomeVersion = 'hg19',
                      outDir = getwd()) {

  if (genomeVersion == 'hg19') {
    species <- 'human'
  } else if (genomeVersion == 'mm9') {
    species <- 'mouse'
  } else if (genomeVersion == 'ce10') {
    species <- 'worm'
    #TODO#
    #msigdb and GO term analysis is not supported for ce10 genome
    #because of gene id issues. However, this will be fixed
    if (msigdbAnalysis == TRUE) {
      warning('Turning off msigdbAnalysis for genomeVersion ce10\n')
      msigdbAnalysis <- FALSE
    }
    if (goAnalysis == TRUE) {
      warning('Turning off msigdbAnalysis for genomeVersion ce10\n')
      goAnalysis <- FALSE
    }
  } else if (genomeVersion == 'dm3') {
    species <- 'fly'
  } else {
    stop (genomeVersion,' is not a supported genome version.')
  }

  if (species != 'human') {
    if (queryFilePath == 'testdata') {
      stop('Test-data only works for human.
           Please provide a queryFilePath to input in BED format \n')
    }

    if (gffFilePath == 'testdata') {
      stop('Test-data only works for human.
           Please provide a gffFilePath to input in GTF format \n')
    }

    if (msigdbFilePath == 'testdata' && msigdbAnalysis == TRUE) {
      stop('Test-data only works for human.
           Please provide a gene set dataset with ENTREZ gene ids
            downloaded from MSIGDB database or
            set msigdbAnalysis option to FALSE \n')
    }
  }

  reportFile <- system.file('report.Rmd', package='RCAS')
  headerFile <- system.file('header.html', package='RCAS')
  outFile <- paste0(basename(queryFilePath), '.RCAS.report.html')

  rmarkdown::render(
    input = reportFile,
    output_dir = outDir,
    output_file = outFile,
    output_format = rmarkdown::html_document(
      toc = TRUE,
      toc_float = TRUE,
      theme = 'simplex',
      number_sections = TRUE,
      includes = rmarkdown::includes(in_header = headerFile)
      ),
    params = list(query = queryFilePath,
                  gff = gffFilePath,
                  msigdb = msigdbFilePath,
                  annotationSummary = annotationSummary,
                  goAnalysis = goAnalysis,
                  msigdbAnalysis = msigdbAnalysis,
                  motifAnalysis = motifAnalysis,
                  genomeVersion = genomeVersion,
                  species = species)
    )
}
