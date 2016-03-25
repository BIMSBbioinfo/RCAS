#' createControlRegions
#'
#' Given a GRanges object of query regions, create a background set of peaks that have
#' the same length distribution based on the flanking regions of the peaks based on the desription in motifRG's paper:  "For each peak in each dataset, we
#' first chose one corresponding background sequence from the flanking regions, randomly
#' chosen from either side <200 nt from the edge of the peak, and with the same width as
#' the peak."
#'
#' @param queryRegions GRanges object containing coordinates of input query
#' regions imported by the \code{\link{importBed}} function.
#' @return GRanges object that contains the same number of regions as query regions
#' @examples
#' bed <- importBed('input.bed')
#' controlRegions <- createControlRegions(queryRegions=bed)
#'
#' @export
createControlRegions <- function (queryRegions) {

  upFlank <- GenomicRanges::flank(queryRegions, start=TRUE, width=200)
  downFlank <- GenomicRanges::flank(queryRegions, start=FALSE, width=200)

  s <- sample(c(0:1), length(queryRegions), replace = TRUE) #randomly select up or downstream flanking regions
  flanks <- upFlank
  flanks[s==1] <- downFlank[s==1]

  df <- data.frame(seqnames = GenomicRanges::seqnames(queryRegions),
                   start    = GenomicRanges::start(queryRegions),
                   end      = GenomicRanges::end(queryRegions),
                   strand   = GenomicRanges::strand(queryRegions),
                   width    = GenomicRanges::width(queryRegions))

  df$flankStart <- GenomicRanges::start(flanks)
  df$flankEnd <- GenomicRanges::end(flanks)
  df$randomFlankStart <- apply(df[,6:7], 1, FUN = function(x) sample(x[1]:x[2], 1))
  df$randomFlankEnd <- df$randomFlankStart + df$width - 1

  controlRegions <- df[1] #control intervals
  controlRegions$start <- df$randomFlankStart
  controlRegions$end <- df$randomFlankEnd
  controlRegions$id <- paste(controlRegions$seqnames, controlRegions$start, controlRegions$end, df$strand, sep=':')
  controlRegions$score <- 1
  controlRegions$strand <- df$strand

  #convert to Granges object
  controlRegions <- makeGRangesFromDataFrame(controlRegions)

  return(controlRegions)
}

#' Extract genomice sequences from BSgenome libraries
#'
#' Given a GRanges object and a genome version (hg19, mm9, ce6, or mm9), this function will
#' extract the DNA sequences for all the genomic regions found in the input object.
#'
#' @param queryRegions GRanges object containing coordinates of input query
#' regions imported by the \code{\link{importBed}} function.
#' @param genomeVersion A character string to denote which BS genome library should be used to
#' extract the sequences. Available options are hg19, mm9, ce6, and mm9
#'
#' @return DNAStringSet object will be returned
#'
#' @examples
#'
#' bed <- importBed('input.bed')
#' sequences <- extractSequences(queryRegions=bed, genomeVersion='hg19')
#'
#' @export
extractSequences <- function (queryRegions, genomeVersion) {

  if (genomeVersion == 'hg19') {
    seqDb = BSgenome.Hsapiens.UCSC.hg19::Hsapiens
  } else if(genome_version == 'mm9') {
    seqDb = BSgenome.Mmusculus.UCSC.mm9::Mmusculus
  } else if(genome_version == 'ce6') {
    seqDb = BSgenome.Celegans.UCSC.ce6::Celegans
  } else if(genome_version == 'dm3') {
    seqDb = BSgenome.Dmelanogaster.UCSC.dm3::Dmelanogaster
  } else {
    stop ("Cannot extract fasta sequences from genome versions except: hg19, mm9, ce6, and dm3\n")
  }
  GenomeInfoDb::seqlevelsStyle(seqDb) =  GenomeInfoDb::seqlevelsStyle(queryRegions)
  sequences <- Biostrings::getSeq(seqDb, queryRegions)
  names(sequences) <- paste('seq',1:length(sequences), sep='_')

  return (sequences)
}

#' runMotifRG
#'
#' This function makes use of \code{motifRG} library to carry out de novo motif discovery
#' from the input query regions
#'
#' @param queryRegions GRanges object containing coordinates of input query
#' regions imported by the \code{\link{importBed}} function.
#' @param genomeVersion A character string to denote which BS genome library should be used to
#' extract the sequences. Available options are hg19, mm9, ce6, and mm9
#' @param motifN A positive integer (default:5) denoting the maximum number of motifs that
#' should be sought by the \code{motifRG::findMotifFgBg} function
#'
#' @return a list of objects returned by the \code{motifRG::findMotif} function
#' @examples
#'
#' bed <- importBed('input.bed')
#' motifResults <- runMotifRG(queryRegions=bed, genomeVersion='hg19')

#' @export
runMotifRG <- function (queryRegions, genomeVersion, motifN = 5) {

  controlRegions <- createControlRegions(queryRegions)

  cat('extracting peak sequences from fasta..\n')
  querySeqs <- extractSequences(queryRegions, genomeVersion)

  cat('extracting background sequences from fasta..\n')
  controlSeqs <- extractSequences(controlRegions, genomeVersion)

  cat('running motifRG...')
  motifResults <- motifRG::findMotifFgBg(querySeqs, controlSeqs, enriched.only=T, max.motif=motifN , both.strand=FALSE)

  return(motifResults)
}

#' getMotifSummaryTable
#'
#' A repurposed/simplified version of the \code{motifRG::summaryMotif} function.
#'
#' @param motifResults Output object of \code{runMotifRG} function
#' @return A data.frame object containing summary statistics about the discovered motifs
#' @examples
#'
#' bed <- importBed('input.bed')
#' motifResults <- runMotifRG(queryRegions=bed, genomeVersion='hg19')
#' motifSummary <- getMotifSummaryTable(motifResults)
#'
#' @export
getMotifSummaryTable <- function(motifResults){

  if(is.null(motifResults)) {return(NULL)}
  motifs <- motifResults$motifs
  category <- motifResults$category
  scores <- c()
  motifPatterns <- c()
  hitsCounts1 <- c()
  hitsCounts2 <- c()
  seqCounts1 <- c()
  seqCounts2 <- c()
  fgSet <- category == 1
  bgSet <- !fgSet
  fgSize <- sum(fgSet)
  bgSize <- sum(bgSet)
  #summarize motif
  for(i in 1:length(motifs)){
    motifPatterns <- c(motifPatterns, motifs[[i]]@pattern)
    scores[i] <- motifs[[i]]@score
    count <- motifs[[i]]@count
    hitsCounts1[i] <- sum(count[fgSet])
    hitsCounts2[i] <- sum(count[bgSet])
    seqCounts1[i] <-  sum(count[fgSet] > 0)
    seqCounts2[i] <-  sum(count[bgSet] > 0)
  }
  ratio <- (hitsCounts1/hitsCounts2)/(fgSize/bgSize)
  frac1 <- seqCounts1/fgSize
  frac2 <- seqCounts2/bgSize
  summary <- data.frame(patterns=motifPatterns,
                        scores=round(scores,1),
                        fgHits=hitsCounts1, bgHits= hitsCounts2,
                        fgSeq=seqCounts1, bgSeq = seqCounts2,
                        ratio=round(ratio,1), fgFrac=round(frac1,4), bgFrac=round(frac2,4))
  return(summary)
}
