#preprocess the gtf/gff files to 
# 1. create GRanges object and write it to .rds file; 
# 2. create txdb file using GenomicFeatures library and write to disk 
# Author: BU 

#1. Collect arguments
args <- commandArgs(TRUE)

## Default setting when no arguments passed
if(length(args) < 1) {
  args <- c("--help")
}

help_command = "
preprocess.anot.R: given a gtf/gff3 file; create rtracklayer::GRanges and GenomicFeatures::txdb objects and write to disk

Arguments:
--gff_file=<path to ENSEMBL GTF/GFF file>  - e.g /data/akalin/Base/Annotation/GenomeAnnotation/hg19/EnsemblGenes/*.gtf
--peak_file=<path to query regions in BED format> - e.g. Hafner2009.bed
--help              - print this text

Example:
Rscript preprocess.anot.R --gff_file=/data/akalin/Base/Annotation/GenomeAnnotation/hg19/EnsemblGenes/142812_EnsemblFTP_hg19_EnsemblGenes.gtf"
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

if(!("gff_file" %in% argsDF$V1)) {
  cat(help_command, "\n")
  stop("provide the path to gtf/gff file")
}

#load libraries#
suppressMessages(library('tools'))
suppressMessages(library('rtracklayer'))
suppressMessages(library('GenomicFeatures'))

gff_file = argsL$gff_file
peak_file = argsL$peak_file

gff_version = file_ext(gff_file)

if(file.exists(paste0(gff_file, ".granges.rds")))
{
  cat('Reading existing granges.rds object\n')
  gff = readRDS(paste0(gff_file, ".granges.rds"))
}else {
  if (gff_version == 'gtf'){
    cat('importing gtf file\n')
    gff = import.gff(gff_file)
  }else if (gff_version == 'gff3'){
    cat('importing gff3 file\n')
    gff = import.gff3(gff_file)
  }else{
    stop('Only gtf and gff3 file extensions are allowed\n')
  } 
  #convert seqlevelsStyle to 'UCSC" 
  seqlevelsStyle(gff) = 'UCSC'
  gff = keepStandardChromosomes(gff)
  saveRDS(gff, file=paste0(gff_file, ".granges.rds"))
}

if(!is.null(gff)){
  if(!file.exists(paste0(gff_file, 'gfeatures.txdb')))
  {  
    cat('Creating txdb object from GRanges\n')
    txdb = makeTxDbFromGRanges(gr = gff)
    cat('Saving txdb into file\n')
    saveDb(txdb, file = paste0(gff_file,'gfeatures.txdb'))
  }else{
    cat('Txdb already exists.. Loading from ',paste0(gff_file, 'gfeatures.txdb'),'\n')
  }
}else{
  stop('imported GFF is empty!\n')
}

##Find all gene ids available in GFF and find genes that overlap with the query regions in bed file
peaks = import.bed(peak_file)
seqlevelsStyle(peaks) = 'UCSC'

all_gene_ids = na.omit(unique(gff$gene_id))
overlaps = gff[queryHits(findOverlaps(gff, peaks))]
targeted_gene_ids = na.omit(unique(overlaps$gene_id))

out_prefix = sub(x = basename(peak_file), pattern = file_ext(peak_file), replacement = '')

write(x = all_gene_ids, file = paste0(out_prefix, '.background_genes.txt'))
write(x = targeted_gene_ids, file = paste0(out_prefix, '.targeted_genes.txt'))




