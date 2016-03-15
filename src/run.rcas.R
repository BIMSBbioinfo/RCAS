#preprocess the gtf/gff files to 
# 1. create GRanges object and write it to .rds file; 
# 2. create txdb file using GenomicFeatures library and write to disk 
# Author: BU 

#1. Collect arguments
args <- commandArgs(TRUE)

## Default setting when no arguments passed
if(length(args) < 3) {
  args <- c("--help")
}

help_command = "
run.rcas.R: run RCAS pipeline for a given ENSEMBL gtf file and BED formatted input file 

Arguments:
--gff_file=<path to ENSEMBL GTF/GFF file>  - e.g /data/akalin/Base/Annotation/hg19/EnsemblGenes/142812_EnsemblFTP_hg19_EnsemblGenes.gtf
--peak_file=<path to query regions in BED format> - e.g. Hafner2009.bed
--genome_version=<choose genome version: e.g. hg19; dm3; ce6; mm9>
--help              - print this text"

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
  stop("provide the path to gtf/gff file\n")
}

if(!("peak_file" %in% argsDF$V1)) {
  cat(help_command, "\n")
  stop("provide the path to input BED file\n")
}

if(!("genome_version" %in% argsDF$V1)) {
  cat(help_command, "\n")
  stop("provide the genome version: choose between hg19, mm9, ce6 or dm3\n")
}

#load libraries#
suppressMessages(library('tools'))
suppressMessages(library('rtracklayer'))
suppressMessages(library('GenomicFeatures'))

gff_file = argsL$gff_file
peak_file = argsL$peak_file
genome_version = argsL$genome_version

### must be configured by configuration script 
source_dir = '/home/buyar/projects/RCAS/src'
base_dir = '/home/buyar/projects/RCAS/src/base'
Rscript = '/usr/bin/Rscript'
###

if(!genome_version %in% c('hg19', 'ce6', 'mm9', 'dm3')){
  cat(help_command,"\n")
  stop("Error: The supported genome versions are hg19, ce6, mm9 and dm3\n")
}

if(genome_version %in% c('hg19', 'hg38')){
  species = 'human'
}else if(genome_version %in% c('mm9', 'mm10')){
  species = 'mouse'
}else if(genome_version %in% c('ce6', 'ce10')){
  species = 'worm'
}else if(genome_version %in% c('dm3', 'dm6')){
  species = 'fly'
}

out_prefix = sub(x = basename(peak_file), pattern = paste0('.',file_ext(peak_file)), replacement = '')


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
  if(!file.exists(paste0(gff_file, '.gfeatures.txdb')))
  {  
    cat('Creating txdb object from GRanges\n')
    txdb = makeTxDbFromGRanges(gr = gff)
    cat('Saving txdb into file\n')
    saveDb(txdb, file = paste0(gff_file,'.gfeatures.txdb'))
  }else{
    cat('Txdb already exists.. Loading from ',paste0(gff_file, '.gfeatures.txdb'),'\n')
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


background_geneset = paste0(out_prefix, '.background_genes.txt')
write(x = all_gene_ids, file = background_geneset)
targeted_geneset = paste0(out_prefix, '.targeted_genes.txt')
write(x = targeted_gene_ids, file = targeted_geneset)

#use the BED file to run rcas.motif.R module
MOTIF_command = paste0(Rscript,' ',source_dir,'/rcas.motif.R',
                       ' --peak_file=',peak_file,
                       ' --genome_version=',genome_version)

cat(MOTIF_command,'\n')
system(MOTIF_command)

#Use the gene lists to run GO term and msigdb analyses 
GO_command = paste0(Rscript,' ',source_dir,'/rcas.GO.R',
                  ' --background_list=',background_geneset,
                  ' --targeted_list=',targeted_geneset,
                   ' --out_prefix=',out_prefix,
                   ' --species=',species)
cat(GO_command,'\n')
system(GO_command)

#Use the gene lists to run MSIGDB module
msigdb_gmt = paste0(base_dir,'/c2.cp.v5.0.entrez.',genome_version,'.gmt')
MSIGDB_command = paste0(Rscript,' ',source_dir,'/rcas.msigdb.R',
                    ' --gmt=',msigdb_gmt,
                    ' --background_list=',background_geneset,
                    ' --targeted_list=',targeted_geneset,
                    ' --out_prefix=',out_prefix,
                    ' --species=',species)
cat(MSIGDB_command,'\n')
system(MSIGDB_command)

work_dir = getwd()
report_script = paste0(source_dir,'/rcas.Rmd')
output_filename = paste0(out_prefix, '.rcas.html')
css = paste0(base_dir,'/custom.css')
header=paste0(base_dir,'/header.html')
REPORT_command = paste0(Rscript," -e \"library('rmarkdown'); rmarkdown::render('",report_script,"',",
                          " output_file = '",output_filename,"',",
                          " intermediates_dir = '",work_dir,"',",
                          " output_dir = '",work_dir,"',",
                          " html_document(toc=TRUE, theme='cerulean', number_sections=TRUE, css='",css,"',",
                          " includes=includes(before_body='",header,"')))\"",
                          " ",work_dir," ",peak_file," ",paste0(gff_file,'.granges.rds')," ",
                          " ",paste0(gff_file,'.gfeatures.txdb'),
                          " ",paste0(out_prefix,".BP.GO.results.tsv"),
                          " ",paste0(out_prefix,".MF.GO.results.tsv"),
                          " ",paste0(out_prefix,".CC.GO.results.tsv"),
                          " ",paste0(out_prefix,".msigdb_results.tsv"),
                          " RUN")

cat(REPORT_command)
system(REPORT_command)
