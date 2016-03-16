# Script to run RCAS in parallel 

#1. Collect arguments
args <- commandArgs(TRUE)

## Default setting when no arguments passed
if(length(args) < 3) {
  args <- c("--help")
}

help_command = "

Arguments:
--gff_file=<path to ENSEMBL GTF/GFF file>  - e.g /data/akalin/Base/Annotation/GenomeAnnotation/hg19/EnsemblGenes/*.gtf
--bed_dir=<path to directory containing input bed files> e.g. /data/akalin/buyar/rcas_analysis/sample_bed_files/
--genome_version=<choose genome version: e.g. hg19; dm3; ce6; mm9>
--help              - print this text

Example:
Rscript run_rcas.parallel.R --gff_file=/data/akalin/Base/Annotation/hg19/EnsemblGenes/142812_EnsemblFTP_hg19_EnsemblGenes.gtf --genome_version=hg19 --bed_dir=./bed_files/"

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
  stop("provide the path to gtf file")
}
if(!("genome_version" %in% argsDF$V1)) {
  cat(help_command, "\n")
  stop("Choose which species' annotation to retrieve: hg19, mm9, dm3, ce6")
}

if(!("bed_dir" %in% argsDF$V1)) {
  cat(help_command, "\n")
  stop("Provide path to directory containing input bed files")
}

## load necessary libraries
library(tools)
library(foreach)
library(doParallel)

gff_file = argsL$gff_file
genome_version = argsL$genome_version
bed_dir = argsL$bed_dir


### must be configured by configuration script 
source_dir = '/home/buyar/projects/RCAS/src'
Rscript = '/usr/bin/Rscript'
###


bedfiles = list.files(path=bed_dir, pattern='*.bed', full.names = TRUE)
workdir = getwd()

##register clusters for parallel executoin
registerDoParallel(cores = 5)

foreach (i=1:length(bedfiles)) %dopar% { 
  setwd(workdir)
  f = bedfiles[i]
  f = file_path_as_absolute(f)
  base = sub(pattern = '.bed$', basename(f), replacement = '')
  cat(f, base, '\n')
  system(paste('mkdir', base))
  setwd(base) 
  rcas_command = paste0(Rscript," ",source_dir,"/run.rcas.R --gff_file=",gff_file," --genome_version=",genome_version," --peak_file=",f)
  cat('workdir=',getwd(),'\n')
  cat('rcascommand',rcas_command,'\n')
  system(rcas_command)
}











