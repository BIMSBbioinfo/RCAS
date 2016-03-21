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
          message('File ',rdsFilePath,'already exists.\n Use overwriteObjectAsRds = TRUE to overwrite the file\n')
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

#' @export
getTxdbFeatures <- function (gff) {

  txdb = GenomicFeatures::makeTxDbFromGRanges(gff)

  transcripts = GenomicFeatures::transcripts(txdb) #has tx_name
  transcripts$gene_name = gff[match(transcripts$tx_name, gff$transcript_id)]$gene_name

  exons = unlist(GenomicFeatures::exonsBy(txdb, by="tx", use.names=TRUE)) #row names are transcript_ids
  exons$tx_name = names(exons)
  exons$gene_name = gff[match(names(exons), gff$transcript_id)]$gene_name

  introns = unlist(GenomicFeatures::intronsByTranscript(txdb, use.names=TRUE))#row names are transcript_ids
  introns$tx_name = names(introns)
  introns$gene_name = gff[match(names(introns), gff$transcript_id)]$gene_name

  promoters = GenomicFeatures::promoters(txdb) #has tx_name
  promoters$gene_name = gff[match(promoters$tx_name, gff$transcript_id)]$gene_name

  fiveUTRs = unlist(GenomicFeatures::fiveUTRsByTranscript(txdb,use.names=TRUE))#row names are transcript_ids
  fiveUTRs$tx_name = names(fiveUTRs)
  fiveUTRs$gene_name = gff[match(names(fiveUTRs), gff$transcript_id)]$gene_name

  threeUTRs = unlist(GenomicFeatures::threeUTRsByTranscript(txdb,use.names=TRUE))#row names are transcript_ids
  threeUTRs$tx_name = names(threeUTRs)
  threeUTRs$gene_name = gff[match(names(threeUTRs), gff$transcript_id)]$gene_name

  cds = unlist(GenomicFeatures::cdsBy(txdb, by="tx",use.names=TRUE))#row names are transcript_ids
  cds$tx_name = names(cds)
  cds$gene_name = gff[match(names(cds), gff$transcript_id)]$gene_name

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
  overlapsQuery = as.data.frame(queryRegions[GenomicRanges::queryHits(overlaps)])
  overlapsGff = gff[GenomicRanges::subjectHits(overlaps)]
  overlapsGff$overlappingQuery = paste(overlapsQuery$seqnames,
                                       overlapsQuery$start,
                                       overlapsQuery$end,
                                       overlapsQuery$strand,
                                       sep=':')
  return (overlapsGff)
}







