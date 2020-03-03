#' createControlRegions
#'
#' Given a GRanges object of query regions, create a background set of peaks
#' that have the same length distribution based on the flanking regions of the
#' peaks.
#'
#' @param queryRegions GRanges object containing coordinates of input query
#'   regions imported by the \code{\link{importBed}} function.
#' @return GRanges object that contains the same number of regions as query
#'   regions
#' @examples
#' data(queryRegions)
#' controlRegions <- createControlRegions(queryRegions = queryRegions)
#'
#' @import GenomicRanges
#' @export
createControlRegions <- function (queryRegions) {

  upFlank <- GenomicRanges::flank(queryRegions, start = TRUE, width=200)
  downFlank <- GenomicRanges::flank(queryRegions, start = FALSE, width=200)

  #randomly select up or downstream flanking regions
  s <- sample(c(0:1), length(queryRegions), replace = TRUE)
  flanks <- upFlank
  flanks[s==1] <- downFlank[s==1]

  df <- data.frame(seqnames = GenomicRanges::seqnames(queryRegions),
                   start    = GenomicRanges::start(queryRegions),
                   end      = GenomicRanges::end(queryRegions),
                   strand   = GenomicRanges::strand(queryRegions),
                   width    = GenomicRanges::width(queryRegions))

  df$flankStart <- GenomicRanges::start(flanks)
  df$flankEnd <- GenomicRanges::end(flanks)
  df$randomFlankStart <- apply(X = df[,6:7], 1,
                               FUN = function(x) sample(x[1]:x[2], 1))
  df$randomFlankEnd <- df$randomFlankStart + df$width - 1

  controlRegions <- df[1] #control intervals
  controlRegions$start <- df$randomFlankStart
  controlRegions$end <- df$randomFlankEnd
  controlRegions$id <- paste(controlRegions$seqnames,
                             controlRegions$start,
                             controlRegions$end,
                             df$strand, sep=':')
  controlRegions$score <- 1
  controlRegions$strand <- df$strand

  #convert to Granges object
  controlRegions <- GenomicRanges::makeGRangesFromDataFrame(controlRegions)

  return(controlRegions)
}

#' extractSequences
#'
#' Given a GRanges object and a genome version (hg19, mm9, ce10 or dm3), this
#' function extracts the DNA sequences for all genomic regions found in an input
#' object.
#'
#' @param queryRegions GRanges object containing coordinates of input query
#'   regions imported by the \code{\link{importBed}} function
#' @param genomeVersion A character string to denote the BS genome library
#'   required to extract sequences. Available options are hg19, mm9, ce10 and
#'   dm3.
#'
#' @return DNAStringSet object will be returned
#'
#' @examples
#'
#' data(queryRegions)
#' sequences <- extractSequences(queryRegions = queryRegions,
#'                              genomeVersion = 'hg19')
#' @importFrom BSgenome.Hsapiens.UCSC.hg19 Hsapiens
#' @importFrom GenomeInfoDb seqlevelsStyle
#' @importFrom Biostrings getSeq
#' @export
extractSequences <- function (queryRegions, genomeVersion) {

  if (genomeVersion == 'hg38') {
    seqDb <- BSgenome.Hsapiens.UCSC.hg38::Hsapiens
  } else if (genomeVersion == 'hg19') {
    seqDb <- BSgenome.Hsapiens.UCSC.hg19::Hsapiens
  } else if (genomeVersion == 'mm10') {
    seqDb <- BSgenome.Mmusculus.UCSC.mm10::Mmusculus
  } else if (genomeVersion == 'mm9') {
    seqDb <- BSgenome.Mmusculus.UCSC.mm9::Mmusculus
  } else if (genomeVersion == 'ce10') {
    seqDb <- BSgenome.Celegans.UCSC.ce10::Celegans
  } else if (genomeVersion == 'dm3') {
    seqDb <- BSgenome.Dmelanogaster.UCSC.dm3::Dmelanogaster
  } else {
    stop ("Cannot extract fasta sequences from genome versions except:
          hg38, hg19, mm10, mm9, ce10 and dm3\n")
  }
  style <- GenomeInfoDb::seqlevelsStyle(queryRegions)
  GenomeInfoDb::seqlevelsStyle(seqDb) <- style
  sequences <- Biostrings::getSeq(seqDb, queryRegions)
  names(sequences) <- paste('seq', 1:length(sequences), sep='_')

  return (sequences)
}


#' Generate K-mers
#' 
#' Given a list of characters, generates all possible fixed length strings 
#' 
#' @param k The length of the strings to be generated 
#' @param letters A character vector 
#' @return Vector of strings 
#' @examples
#' generateKmers(3, c('A', 'C', 'G'))
#' 
#' @export
generateKmers <- function(k, letters = c("A", "C", "G", "T")) {
  kmer <- c()
  for(i in 1:k){
    kmer <- unlist(lapply(letters, function(x){paste(kmer, x, sep="")}))
  }
  return(kmer)
}

countPattern <- function(seqs, patterns, maxMismatch = 0, nCores = 1) {
  cl <- parallel::makeCluster(nCores)
  parallel::clusterExport(cl = cl, varlist = c('seqs', 'patterns', 'maxMismatch'), 
                          envir = environment())
  M <- do.call(rbind, pbapply::pblapply(cl = cl, X = patterns, 
                               FUN = function(x) {
                                 Biostrings::vcountPattern(x, seqs, 
                                                           max.mismatch = maxMismatch)
                               }))
  parallel::stopCluster(cl)
  rownames(M) <- patterns
  return(t(M))
}

extractMatches <- function(seqs, patterns, minMismatch = 0, maxMismatch = 0) {
  matches <- lapply(patterns, function(x) {
    m <- Biostrings::vmatchPattern(pattern = x, subject = seqs, max.mismatch = maxMismatch)
    m <- unlist(m[lengths(m) > 0])
    if(length(m) == 0) {
      return(NULL)
    }
    # remove out of bounds matches
    m <- m[start(m) >= 1]
    m <- m[end(m) <= IRanges::width(seqs[names(m)])]
    m <- split(m, names(m))
    common <- intersect(names(m), names(seqs))
    if(length(common) == 0) {
      return(NULL)
    }
    paste(unlist(Biostrings::extractAt(seqs[common], m[common])))
  })
  names(matches) <- patterns
  return(matches)
}

#' Find Differential Motifs 
#' 
#' @export
findDifferentialMotifs <- function(querySeqs, 
                                   controlSeqs, 
                                   motifWidth = 6,
                                   motifN, 
                                   nCores, 
                                   maxMismatch = 1) {
  
  kmers <- generateKmers(width = motifWidth)
  
  selected <- seq(1:length(querySeqs))
  if(length(querySeqs) > 1000) {
    selected <- sample(1:length(querySeqs), 1000)
  }
  
  query <- countPattern(querySeqs[selected], kmers, maxMismatch, nCores)
  ctrl <- countPattern(controlSeqs[selected], kmers, maxMismatch, nCores)

  # only keep patterns that are more frequent in the query
  queryHits <- apply(query, 2, function(x) sum(x > 0))
  controlHits <- apply(ctrl, 2, function(x) sum(x > 0))
  candidates <- names(which(log2((queryHits + 1) / (controlHits+1)) > 0))
  
  query <- query[,candidates]
  ctrl <- ctrl[,candidates]
  
  # train a random forest model 
  df <- as.data.frame(rbind(query[,candidates], ctrl[,candidates]))
  df$label <- c(rep("query", nrow(query)), 
                   rep("ctrl", nrow(ctrl)))
  
  fit <- ranger::ranger(label ~ ., df, importance = 'impurity')
  var.imp <- sort(ranger::importance(fit), decreasing = T)
  #get top variables 
  top <- names(var.imp[1:motifN])

  results <- list("counts_query" = query[,top, drop = F], 
                  "counts_ctrl" = ctrl[,top,drop = F], 
                  "matches_query" = extractMatches(seqs = querySeqs, 
                                                   patterns = top, 
                                                   maxMismatch = maxMismatch), 
                  "matches_ctrl" = extractMatches(seqs = controlSeqs, 
                                                  patterns = top, 
                                                  maxMismatch = maxMismatch))
  return(results)
}

#' run motifRG 
#' 
#' @export
runMotifRG <- function() {
  .Deprecated("runMotifDiscovery")
}

#' runMotifDiscovery
#'
#' This function builds a random forest classifier to find the top most
#' discriminative motifs in the query regions compared to the background. The
#' background sequences are automatically generated based on the query regions.
#' First k-mers of a fixed length are generated. The query and control sequences
#' are searched for k-mers allowing for mismatches. A random forest model is 
#' trained to find the most discriminative motifs.  
#' 
#' @param queryRegions GRanges object containing coordinates of input query 
#'   regions imported by the \code{\link{importBed}} function
#' @param resizeN Integer value (default: 0) to resize query regions if they are
#'   shorter than the value of \code{resize}. Set to 0 to disable resize.
#' @param sampleN A positive integer value. The queryRegions are randomly 
#'   downsampled to include intervals as many as \code{sampleN}. The input will
#'   be downsampled only if this value is larger than zero and less than the
#'   total number of input intervals.
#' @param genomeVersion A character string to denote the BS genome library 
#'   required to extract sequences. Available options are hg19, mm9, ce10 and 
#'   dm3.
#' @param motifN A positive integer (default:5) denoting the maximum number of 
#'   motifs that should be sought by the \code{motifRG::findMotifFgBg} function
#' @param nCores A positive integer (default:4) number of cores used for 
#'   parallel execution.
#' @return a list of objects returned by the \code{motifRG::findMotif} function
#' @examples
#' data(queryRegions)
#' motifResults <- runMotifDiscovery(queryRegions = queryRegions[1:1000],
#'                                   genomeVersion = 'hg19',
#'                                   motifWidth = 6, 
#'                                   resize = 15,
#'                                   motifN = 1,
#'                                   maxMismatch = 1,
#'                                   nCores = 2)
#' @import IRanges
#' @export
runMotifDiscovery <- function (queryRegions, resizeN = 0, motifWidth = 6,
                        sampleN = 0, genomeVersion, maxMismatch = 1,
                        motifN = 5, nCores = 4) {

  if(sampleN > 0 && length(queryRegions) > sampleN) {
    message("Randomly sampling query regions for motif analysis. 
            Downsampling to ",sampleN," regions")
    queryRegions <- sample(queryRegions, sampleN)
  } 
  
  if(resizeN > 0) {
    resizeIntervals <- width(queryRegions) < resizeN
    if(sum(resizeIntervals) > 0){
      message("Found ",sum(resizeIntervals)," query regions shorter than ", 
              resizeN, " bps. Resizing those regions to ",resizeN,' bps')
      queryRegions[resizeIntervals] <- GenomicRanges::resize(
        x = queryRegions[resizeIntervals], 
        width = 15, 
        fix = 'center')
    }
  }

  message('extracting sequences from fasta..')
  querySeqs <- extractSequences(queryRegions, genomeVersion)

  message('extracting background sequences from fasta..')
  controlRegions <- createControlRegions(queryRegions)
  controlSeqs <- extractSequences(controlRegions, genomeVersion)

  message('running motif discovery ... ')
  motifResults <- findDifferentialMotifs(querySeqs = querySeqs, 
                                         controlSeqs = controlSeqs, 
                                         motifWidth = 6, 
                                         motifN = motifN,
                                         maxMismatch = maxMismatch, 
                                         nCores = nCores)
  return(motifResults)
}

#' getMotifSummaryTable
#'
#' Get summary stats for top discovered motifs 
#'
#' @param motifResults Output object of \code{runMotifDiscovery} function
#' @return A data.frame object containing summary statistics about the
#'   discovered motifs
#' @examples
#'
#' data(queryRegions)
#' motifResults <- runMotifDiscovery(queryRegions = queryRegions[1:1000],
#'                                   genomeVersion = 'hg19',
#'                                   resize = 15,
#'                                   motifN = 1,
#'                                   maxMismatch = 1,
#'                                   nCores = 2)
#' motifSummary <- getMotifSummaryTable(motifResults)
#'
#' @export
getMotifSummaryTable <- function(motifResults){

  data.frame('patterns' = names(motifResults$matches_query),
             'queryHits' = colSums(motifResults$counts_query),
             'controlHits' = colSums(motifResults$counts_ctrl),
             'querySeqs' = colSums(motifResults$counts_query > 0),
             'controlSeqs' = colSums(motifResults$counts_ctrl > 0),
             'queryFraction' = round(colSums(motifResults$counts_query > 0) / 
                                       nrow(motifResults$counts_query), 2),
             'controlFraction' = round(colSums(motifResults$counts_ctrl > 0) / 
                                         nrow(motifResults$counts_ctrl), 2))
}


#' discoverFeatureSpecificMotifs
#' 
#' This function groups query regions based on their overlap with different
#' transcript features and generates a table of top enriched motif and matching
#' patterns for each given transcript feature type along with some other motif 
#' discovery related statistics.
#' 
#' @param queryRegions GRanges object containing coordinates of input query
#'   regions imported by the \code{\link{importBed}} function
#' @param txdbFeatures A list of GRanges objects where each GRanges object
#'   corresponds to the genomic coordinates of gene features such as promoters,
#'   introns, exons, 5'/3' UTRs and whole transcripts.
#'   This list of GRanges objects are obtained by the function
#'   \code{\link{getTxdbFeaturesFromGRanges}} or \code{\link{getTxdbFeatures}}.
#' @param ... Other arguments passed to \code{\link{runMotifRG}} function.
#'   Important arguments are 'genomeVersion' and motifN. If motifN is bigger
#'   than 1, then multiple motifs will be found but only the top motif will be
#'   plotted.
#' @examples
#' \dontrun{
#' data(gff)
#' data(queryRegions)
#' txdbFeatures <- getTxdbFeaturesFromGRanges(gffData = gff)
#' discoverFeatureSpecificMotifs(queryRegions = queryRegions, 
#' genomeVersion = 'hg19', txdbFeatures = txdbFeatures, 
#' motifN = 1, nCores = 1)}
#' 
#' @return A data.frame object  
#' @export 
discoverFeatureSpecificMotifs <- function(queryRegions, txdbFeatures, ...) {
  results <- do.call(rbind, lapply(names(txdbFeatures), function(f) {
    message("Looking for motifs in feature:",f)
    featureCoords <- txdbFeatures[[f]]
    #find query regions that overlap the target features
    q <- queryRegions[unique(queryHits(findOverlaps(queryRegions, featureCoords)))]
    
    #notice that this must be in line with the select statement
    #from getMotifSummaryTable function's output
    motifStats <- data.frame('patterns' = 'NONE', 
                             'queryHits' = 0, 
                             'controlHits' = 0,
                             'querySeqs' = 0,
                             'controlSeqs' = 0,
                             'queryFraction' = 0,
                             'controlFraction' = 0,
                             'feature' = f)
    
    if(length(q) > 0) {
      motifResults <- runMotifDiscovery(queryRegions = q, ...)
      motifStats <- subset(getMotifSummaryTable(motifResults), 
                           select = c('patterns', 
                                      'queryFraction', 'controlFraction'))
      motifStats$querySeqTotal <- length(q)
      motifStats$feature <- f
      }
    return(motifStats)
  }))  
  return(results)
}
