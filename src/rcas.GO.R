#1. Collect arguments
args <- commandArgs(TRUE)

## Default setting when no arguments passed
if(length(args) < 2) {
  args <- c("--help")
}

help_command = "
      rcas.GO.R: GO Enrichment Analysis module of RCAS 

      Arguments:
      --gff3=<path to GENCODE Annotation File>  - e.g /data/akalin/Base/Annotation/GenomeAnnotation/hg19/gencode/gencode.v19.annotation.gff3
      --anot=<path to parse_anot.py output>  - e.g /home/buyar/projects/RCAS/test/PARCLIP_AGO1234_Hafner2010a_hg19_xaa.anot.tsv
      --help              - print this text
      
      Example:
      Rscript rcas.GO.R --gff3=gencode.v19.annotation.gff3 --anot=PARCLIP_AGO1234_Hafner2010a_hg19_xaa.anot.tsv"

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

## Arg1 default
if(!("gff3" %in% argsDF$V1)) {
  cat(help_command, "\n")
  stop("provide the path to gff3 file")
}

## Arg2 default
if(!("anot" %in% argsDF$V1)) {
  cat(help_command, "\n")
  stop("provide the path to anot.bed file")
}

gff3_file = argsL$gff3
anot_file = argsL$anot

########################################################################################################################################
suppressWarnings(suppressMessages(library('rtracklayer')))
#2. Get the list of all protein-coding genes from the GTF file 

if(file.exists(paste0(gff3_file, ".rds")))
{
  gff = readRDS(paste0(gff3_file, ".rds"))
}else
{
  gff = import.gff3(gff3_file)
  saveRDS(gff, file=paste0(gff3_file, ".rds"))
}
idx <- mcols(gff)$gene_type == "protein_coding" & mcols(gff)$type=="gene"
all_gene_ids = unique(gff[idx]$gene_id)
head(all_gene_ids)

#3. Get the list of protein-coding genes from the anot.bed file => those that are found to be targeted in the PAR-CLIP experiments
suppressWarnings(suppressMessages(library(data.table)))

ann = fread(anot_file)
ann.gene_ids = unique(ann$gene_id)

########################################################################################################################################
#4. Look for enriched GO terms for the genes found in the anot.bed compared to the background genes (found in gtf)
#source("http://www.bioconductor.org/biocLite.R")
#biocLite(c("biomaRt", "topGO", "org.Hs.eg.db"))
suppressWarnings(suppressMessages(library('biomaRt')))
suppressWarnings(suppressMessages(library('org.Hs.eg.db'))) #mouse: org.Mm.eg.db
suppressWarnings(suppressMessages(library('topGO')))    

gene_universe = data.table(cbind(all_gene_ids))
colnames(gene_universe) = c('gene_id')
gene_universe$targeted = 0

gene_universe[gene_id %in% ann.gene_ids]$targeted = 1

all_genes = gene_universe$targeted
names(all_genes) = gsub("\\..*$", "", gene_universe$gene_id) #ensembl gene ids don't have trailing version numbers of gene ids: e.g. (ENSG00121311.9 is converted to ENSG00121311) 

selection_function = function (x){ return(x == 1) }

GOdata <- new("topGOdata", ontology = "BP", allGenes = all_genes, geneSel = selection_function, description = "Test", annot = annFUN.org, mapping="org.Hs.eg.db", ID="Ensembl")

saveRDS(GOdata, file="godata.rds")

resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
mt = GenTable(GOdata, classicFisher = resultFisher, topNodes = 10)









