createControlRegions <- function (queryRegions) {

  upFlank <- GenomicRanges::flank(queryRegions, start=TRUE, width=200)
  downFlank <- GenomicRanges::flank(queryRegions, start=FALSE, width=200)

  s <- sample(c(0:1), length(queryRegions), replace = TRUE) #randomly select up or downstream flanking regions
  flanks <- upFlank
  flanks[s==1] <- downFlank[s==1]

  df <- data.frame(chr <- GenomicRanges::seqnames(queryRegions),
                   start <- GenomicRanges::start(queryRegions),
                   end <- GenomicRanges::end(queryRegions),
                   strand <- GenomicRanges::strand(queryRegions),
                   width <- GenomicRanges::width(queryRegions))

  df$flankStart <- GenomicRanges::start(flanks)
  df$flankEnd <- GenomicRanges::end(flanks)
  df$randomFlankStart <- apply(df[,6:7], 1, FUN = function(x) sample(x[1]:x[2], 1))
  df$randomFlankEnd <- df$randomFlankStart + df$width - 1

  controlRegions <- df[1] #control intervals
  controlRegions$start <- df$randomFlankStart
  controlRegions$end <- df$randomFlankEnd
  controlRegions$id <- paste(controlRegions$chr, controlRegions$start, controlRegions$end, df$strand, sep=':')
  controlRegions$score <- 1
  controlRegions$strand <- df$strand

  #convert to Granges object
  controlRegions <- makeGRangesFromDataFrame(controlRegions)

  return(controlRegions)
}

runMotifRG <- function (queryRegions, genomeVersion) {

  controlRegions <- createControlRegions(queryRegions)

  cat('extracting peak sequences from fasta..\n')
  peak_seqs <- getSeq(seq_fasta, peaks)
  names(peak_seqs) <- paste('seq',1:length(peak_seqs), sep='_')

  cat('extracting background sequences from fasta..\n')
  control_seqs <- getSeq(seq_fasta, control_peaks)
  names(control_seqs) <- paste('seq',1:length(control_seqs), sep='_')

  cat('running motifRG...')
  motif_results <- findMotifFgBg(peak_seqs, control_seqs, enriched.only=T, max.motif=5, is.parallel=TRUE, mc.cores=8, both.strand=FALSE)
}


