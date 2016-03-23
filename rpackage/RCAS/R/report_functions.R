#' @export
importGtf <- function (filePath, readFromRds = TRUE, saveObjectAsRds = TRUE, overwriteObjectAsRds = FALSE, keepStandardChr = TRUE) {

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
          message('File ',rdsFilePath,' already exists.\n Use overwriteObjectAsRds = TRUE to overwrite the file\n')
        }
      }
    return (gff)
    }
  } else{
    stop("Couldn't import the gtf from",filePath,"\n")
  }
}

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
    stop("Cannot import BED file because it does not exist at the given location:",filePath,"\n")
  }

}

processHits <- function(queryRegions, tx, type) {
  overlaps <- GenomicRanges::findOverlaps(queryRegions, tx)
  overlapsQuery <- queryRegions[GenomicRanges::queryHits(overlaps)]
  overlapsTX <- tx[GenomicRanges::subjectHits(overlaps)]
  overlapsTX$overlappingQuery <- paste(GenomicRanges::seqnames(overlapsQuery),
                                GenomicRanges::start(overlapsQuery),
                                GenomicRanges::end(overlapsQuery),
                                GenomicRanges::strand(overlapsQuery),
                                sep=':')
  dt <- data.table::data.table('tx_name' = overlapsTX$tx_name, 'query' = overlapsTX$overlappingQuery)
  summary <- dt[,length(unique(query)), by=tx_name]
  colnames(summary) <- c('tx_name', type)
  return(summary)
}

#' @export
getTargetedGenesTable <- function (queryRegions, txdbFeatures) {

  tbls <- lapply(X=seq_along(txdbFeatures),
               FUN=function(i) { processHits(queryRegions = queryRegions,
                                 tx = txdbFeatures[[names(txdbFeatures)[i]]],
                                 type = names(txdbFeatures)[i])})

  tbls <- lapply(tbls, function(i) setkey(i, tx_name))
  merged <- Reduce(function(...) merge(..., all = T), tbls)
  merged$gene_name = gff[match(merged$tx_name, gff$transcript_id)]$gene_name
  merged[is.na(merged)] <- 0
  return(merged)
}

#' @export
getTxdbFeatures <- function (txdb) {

  transcripts = GenomicFeatures::transcripts(txdb) #has tx_name
  exons = GenomicRanges::unlist(GenomicFeatures::exonsBy(txdb, by="tx", use.names=TRUE)) #row names are transcript_ids
  introns = GenomicRanges::unlist(GenomicFeatures::intronsByTranscript(txdb, use.names=TRUE))#row names are transcript_ids
  exonIntronBoundaries = c(GenomicRanges::flank(introns[GenomicRanges::width(introns) >= 100], 50, start=TRUE, both=TRUE),
                           GenomicRanges::flank(introns[GenomicRanges::width(introns) >= 100], 50, start=FALSE, both=TRUE))
  promoters = GenomicFeatures::promoters(txdb) #has tx_name
  fiveUTRs = GenomicRanges::unlist(GenomicFeatures::fiveUTRsByTranscript(txdb,use.names=TRUE))#row names are transcript_ids
  threeUTRs = GenomicRanges::unlist(GenomicFeatures::threeUTRsByTranscript(txdb,use.names=TRUE))#row names are transcript_ids
  cds = GenomicRanges::unlist(GenomicFeatures::cdsBy(txdb, by="tx",use.names=TRUE))#row names are transcript_ids

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


#' @export
getTxdbFeaturesFromGff <- function (gff) {

  txdb = GenomicFeatures::makeTxDbFromGRanges(gff)

  transcripts = GenomicFeatures::transcripts(txdb) #has tx_name
  transcripts$gene_name = gff[match(transcripts$tx_name, gff$transcript_id)]$gene_name

  exons = GenomicRanges::unlist(GenomicFeatures::exonsBy(txdb, by="tx", use.names=TRUE)) #row names are transcript_ids
  exons$tx_name = names(exons)
  exons$gene_name = gff[match(names(exons), gff$transcript_id)]$gene_name

  introns = GenomicRanges::unlist(GenomicFeatures::intronsByTranscript(txdb, use.names=TRUE))#row names are transcript_ids
  introns$tx_name = names(introns)
  introns$gene_name = gff[match(names(introns), gff$transcript_id)]$gene_name

  exonIntronBoundaries = c(GenomicRanges::flank(introns[GenomicRanges::width(introns) >= 100], 50, start=TRUE, both=TRUE),
                           GenomicRanges::flank(introns[GenomicRanges::width(introns) >= 100], 50, start=FALSE, both=TRUE))

  promoters = GenomicFeatures::promoters(txdb) #has tx_name
  promoters$gene_name = gff[match(promoters$tx_name, gff$transcript_id)]$gene_name

  fiveUTRs = GenomicRanges::unlist(GenomicFeatures::fiveUTRsByTranscript(txdb,use.names=TRUE))#row names are transcript_ids
  fiveUTRs$tx_name = names(fiveUTRs)
  fiveUTRs$gene_name = gff[match(names(fiveUTRs), gff$transcript_id)]$gene_name

  threeUTRs = GenomicRanges::unlist(GenomicFeatures::threeUTRsByTranscript(txdb,use.names=TRUE))#row names are transcript_ids
  threeUTRs$tx_name = names(threeUTRs)
  threeUTRs$gene_name = gff[match(names(threeUTRs), gff$transcript_id)]$gene_name

  cds = GenomicRanges::unlist(GenomicFeatures::cdsBy(txdb, by="tx",use.names=TRUE))#row names are transcript_ids
  cds$tx_name = names(cds)
  cds$gene_name = gff[match(names(cds), gff$transcript_id)]$gene_name

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

#' @export
summarizeQueryRegions <- function(queryRegions, txdbFeatures) {
  results = lapply(X = txdbFeatures, FUN = function (x) { length(unique(GenomicRanges::queryHits(GenomicRanges::findOverlaps(queryRegions, x)))) } )
  df <- t(data.frame(results))
  colnames(df) <- c('count')
  return(df)
}

#' @export
queryGff <- function(queryRegions, gff) {
  overlaps = GenomicRanges::findOverlaps(queryRegions, gff) #find all overlapping pairs of intervals between the gff features and BED file (peaks)
  overlapsQuery = queryRegions[GenomicRanges::queryHits(overlaps)]
  overlapsGff = gff[GenomicRanges::subjectHits(overlaps)]
  overlapsGff$overlappingQuery = paste(GenomicRanges::seqnames(overlapsQuery),
                                       GenomicRanges::start(overlapsQuery),
                                       GenomicRanges::end(overlapsQuery),
                                       GenomicRanges::strand(overlapsQuery),
                                       sep=':')
  return (overlapsGff)
}

#' @export
calculateCoverageProfile = function (queryRegions, targetRegions, sampleN = 0){
  windows <- targetRegions[GenomicRanges::width(targetRegions) >= 100]#remove windows shorter than 100 bp

  if (length(windows) > 0) {
    if (sampleN > 0 && sampleN < length(windows)) {
      windows <- windows[sample(length(windows), sampleN)]
    }
    sm <- genomation::ScoreMatrixBin(target = queryRegions, windows = windows, bin.num = 100, strand.aware = TRUE)
    mdata <- as.data.frame(colSums(sm))
    mdata$bins <- c(1:100)
    colnames(mdata) <- c('coverage', 'bins')
    return(mdata)
    } else {
    stop("Cannot compute coverage profile for target regions.\nThere are no target regions longer than 100 bp\n")
  }
}

#' @export
calculateCoverageProfileList <- function (queryRegions, targetRegionsList, sampleN = 0) {
  lapply(X = targetRegionsList, FUN=function(x) { calculateCoverageProfile(queryRegions, x, sampleN = sampleN) })
}

#' @export
calculateCoverageProfileListFromTxdb <- function (queryRegions, txdb, sampleN = 0) {
  txdbFeatures = getTxdbFeatures(txdb = txdb)
  lapply(X = txdbFeatures, FUN=function(x) { calculateCoverageProfile(queryRegions, x, sampleN = sampleN) })
}

#' @export
calculateCoverageProfileFromTxdb <- function (queryRegions, txdb, type, sampleN = 0) {

  if (type == 'transcripts') {
    targetRegions <- GenomicFeatures::transcripts(txdb) #has tx_name
  } else if (type == 'exons') {
    targetRegions = GenomicRanges::unlist(GenomicFeatures::exonsBy(txdb, by="tx", use.names=TRUE)) #row names are transcript_ids
  } else if (type == 'introns') {
    targetRegions <- GenomicRanges::unlist(GenomicFeatures::intronsByTranscript(txdb, use.names=TRUE))#row names are transcript_ids
  } else if (type == 'exonIntronBoundaries') {
    introns <- GenomicRanges::unlist(GenomicFeatures::intronsByTranscript(txdb, use.names=TRUE))#row names are transcript_ids
    targetRegions <- c(GenomicRanges::flank(introns[GenomicRanges::width(introns) >= 100], 50, start=TRUE, both=TRUE),
                             GenomicRanges::flank(introns[GenomicRanges::width(introns) >= 100], 50, start=FALSE, both=TRUE))
  } else if (type == 'promoters') {
    targetRegions <- GenomicFeatures::promoters(txdb) #has tx_name
  } else if (type == 'fiveUTRs') {
    targetRegions = GenomicRanges::unlist(GenomicFeatures::fiveUTRsByTranscript(txdb,use.names=TRUE))#row names are transcript_ids
  } else if (type == 'threeUTRs') {
    targetRegions = GenomicRanges::unlist(GenomicFeatures::threeUTRsByTranscript(txdb,use.names=TRUE))#row names are transcript_ids
  } else if (type == 'cds') {
    targetRegions = GenomicRanges::unlist(GenomicFeatures::cdsBy(txdb, by="tx",use.names=TRUE))#row names are transcript_ids
  } else {
    stop ("Can calculate coverage profiles for only: transcripts, exons, introns, exonIntronBoundaries, promoters, fiveUTRs, threeUTRs, and cds")
  }
  result = calculateCoverageProfile(queryRegions = queryRegions, targetRegions = targetRegions, sampleN = sampleN)
  return (result)
}




