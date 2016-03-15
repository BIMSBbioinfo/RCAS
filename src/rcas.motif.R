# Run motifRG and rGADEM on a given set of BED files 

#1. Collect arguments
args <- commandArgs(TRUE)

## Default setting when no arguments passed
if(length(args) < 2) {
  args <- c("--help")
}

help_command = "
RCAS.motif.R: 

Arguments:
--peak_file=<path to query regions in BED format> - e.g. Hafner2009.bed
--genome_version=<choose genome version: e.g. hg19; dm3; ce6; mm9>
--help              - print this text

Example:
Rscript --peak_file Hafner2009.bed --species --genome_version"
## Help section
if("--help" %in% args) {
  cat(help_command, "\n")
  q(save="no")
}

## Parse arguments (we expect the form --arg=value)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
argsDF <- as.data.frame(do.call("rbind", parseArgs(args)))
argsL <- as.list(as.character(argsDF$V2))
names(argsL) <- argsDF$V1

if(!("peak_file" %in% argsDF$V1)) {
  cat(help_command, "\n")
  stop("provide the path to input BED file")
}

if(!("genome_version" %in% argsDF$V1)) {
  cat(help_command, "\n")
  stop("provide genome version")
}


#load libraries#
suppressMessages(library('tools'))
suppressMessages(library('rtracklayer'))
suppressMessages(library('GenomicFeatures'))
suppressMessages(library('motifRG'))
#suppressMessages(library('rGADEM'))


peak_file = argsL$peak_file
genome_version = argsL$genome_version

if(genome_version == 'hg19'){
  require('BSgenome.Hsapiens.UCSC.hg19')
  seq_fasta = Hsapiens
  seq_fasta
  seqlevelsStyle(seq_fasta) = 'UCSC'
}

#given a list of peaks, create a background set of peaks that have the same length distribution based on the flanking regions of the peaks.

#based on the desription in motifRG's paper:  "For each peak in each dataset, we
#first chose one corresponding background sequence from the flanking
#regions, randomly chosen from either side 0???200 nt from the edge of the
#peak, and with the same width as the peak."

create_background_peaks = function(peaks, out_prefix){
  df = data.frame(chr = seqnames(peaks), start=start(peaks), end=end(peaks), strand=strand(peaks), width=width(peaks))
  up_flank = flank(peaks, start=TRUE, width=200)
  down_flank = flank(peaks, start=FALSE, width=200)
  
  s = sample(c(0:1), length(peaks), replace = TRUE) #randomly select up or downstream flanking regions
  flanks = up_flank
  flanks[s==1] = down_flank[s==1]
  
  df$flank_start = start(flanks)
  df$flank_end = end(flanks)
  
  df$random_flank_start = apply(df[,6:7], 1, FUN = function(x) sample(x[1]:x[2], 1))
  df$random_flank_end = df$random_flank_start + df$width - 1
  
  df2 = df[1] #control intervals
  df2$start = df$random_flank_start
  df2$end = df$random_flank_end
  df2$id = paste(df2$chr, df2$start, df2$end, df$strand, sep=':')
  df2$score = 1
  df2$strand = df$strand
  outfile = paste0(out_prefix, '.background.bed')
  write.table(df2, file = outfile, row.names=FALSE, col.names=FALSE, quote = FALSE, sep='\t')
  return(outfile)
}

cat('importing peak data from',peak_file,'\n')
peaks = import.bed(peak_file)
seqlevelsStyle(peaks) = 'UCSC'

out_prefix = sub(pattern = paste0('.', file_ext(peak_file)), x = basename(peak_file), replacement = '')

cat('creating corresponding background peak data...\n')
outfile = create_background_peaks(peaks, out_prefix)
control_peaks = import.bed(outfile)

#get genome sequences of the corresponding query regions 
cat('extracting peak sequences from fasta..\n')
peak_seqs = getSeq(seq_fasta, peaks)
names(peak_seqs) = paste('seq',1:length(peak_seqs), sep='_')


cat('extracting background sequences from fasta..\n')
control_seqs = getSeq(seq_fasta, control_peaks)
names(control_seqs) = paste('seq',1:length(control_seqs), sep='_')


#write fasta file to disk
peak_fasta = paste0(out_prefix, '.foreground.fasta')
writeXStringSet(peak_seqs, filepath=peak_fasta, format='fasta')

control_peak_fasta = paste0(out_prefix, '.background.fasta')
writeXStringSet(peak_seqs, filepath=control_peak_fasta, format='fasta')

#using motifRG library

cat('running motifRG...')
motif_results = findMotifFgBg(peak_seqs, control_seqs, enriched.only=T,max.motif=5, is.parallel=TRUE, mc.cores=8, both.strand=FALSE) 
motifHtmlTable(motif_results, prefix=out_prefix)


 








