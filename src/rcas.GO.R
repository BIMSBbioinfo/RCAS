#1. Collect arguments
args <- commandArgs(TRUE)

## Default setting when no arguments passed
if(length(args) < 4) {
  args <- c("--help")
}

help_command = "
      rcas.GO.R: GO Enrichment Analysis module of RCAS 

      Arguments:
      --background_list=<background list of genes >  - e.g {sample}.background_genes.txt
      --targeted_list=<peaks file in BED format>  - e.g {sample}.targeted_genes.txt
      --out_prefix=<output prefix>   -e.g Hafner2010.hg19
      --species=<species name>    -Choose between (human, fly, worm, mouse)
      --help              - print this text
      
      Example:
      Rscript rcas.GO.R --background_list={sample}.background_genes.txt--targeted_list={sample}.targeted_genes.txt --out_prefix=Hafner2010.hg19 --species=human"

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
if(!("background_list" %in% argsDF$V1)) {
  cat(help_command, "\n")
  stop("provide the path to background gene list")
}

## Arg2 default
if(!("targeted_list" %in% argsDF$V1)) {
  cat(help_command, "\n")
  stop("provide the path to query/targeted gene list")
}

## Arg3 default
if(!("out_prefix" %in% argsDF$V1)) {
  cat(help_command, "\n")
  stop("provide the output prefix")
}

## Arg4
if(!("species" %in% argsDF$V1)) {
  cat(help_command, "\n")
  stop("Choose which species' annotation to retrieve: human, fly, worm, or mouse")
}

background_list = argsL$background_list
targeted_list = argsL$targeted_list
out_prefix = argsL$out_prefix
species = argsL$species

if (!(species %in% c('human', 'fly', 'worm', 'mouse'))){
  cat(help_command, "\n")
  cat('selected species \"',species, '\" is not supported','\n')
  stop("Choose which species' annotation to retrieve. Supported species are: human, fly, worm, or mouse")
}

########################################################################################################################################
suppressWarnings(suppressMessages(library('rtracklayer')))
#2. Get the list of all protein-coding genes from the GTF file 

all_gene_ids = na.omit(scan(background_list, what = character()))
head(all_gene_ids)

#3. Get the list of genes targeted in the PAR-CLIP experiments
suppressWarnings(suppressMessages(library(data.table)))
targeted_gene_ids = na.omit(scan(targeted_list, what = character()))

########################################################################################################################################
#4. Look for enriched GO terms for the genes found in the anot.bed compared to the background genes (found in gtf)
#source("http://www.bioconductor.org/biocLite.R")
#biocLite(c("biomaRt", "topGO", "org.Hs.eg.db", "org.Ce.eg.db", "org.Mm.eg.db", "org.Dm.eg.db"))
suppressWarnings(suppressMessages(library('biomaRt')))
suppressWarnings(suppressMessages(library('org.Hs.eg.db'))) #mouse: org.Mm.eg.db
suppressWarnings(suppressMessages(library('org.Ce.eg.db'))) 
suppressWarnings(suppressMessages(library('org.Mm.eg.db'))) 
suppressWarnings(suppressMessages(library('org.Dm.eg.db'))) 
suppressWarnings(suppressMessages(library('topGO')))    

gene_universe = data.table(cbind(all_gene_ids))
colnames(gene_universe) = c('gene_id')
gene_universe$targeted = 0

gene_universe[gene_id %in% targeted_gene_ids]$targeted = 1

all_genes = gene_universe$targeted
names(all_genes) = gsub("\\..*$", "", gene_universe$gene_id) #ensembl gene ids don't have trailing version numbers of gene ids: e.g. (ENSG00121311.9 is converted to ENSG00121311) 

selection_function = function (x){ return(x == 1) }

calculate_go_enrichment = function (ontology){

if (species == 'human'){
  species_orgdb = 'org.Hs.eg.db'
}else if (species == 'mouse'){
  species_orgdb = 'org.Mm.eg.db'
}else if (species == 'fly'){
  species_orgdb = 'org.Dm.eg.db'
}else if (species == 'worm'){
  species_orgdb = 'org.Ce.eg.db'
}
  
GOdata <- new("topGOdata", ontology = ontology, allGenes = all_genes, geneSel = selection_function, nodeSize=20, description = "Test", annot = annFUN.org, mapping=species_orgdb, ID="Ensembl")

resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
mt = GenTable(GOdata, classicFisher = resultFisher, topNodes = length(usedGO(GOdata)))

#Do multiple-testing correction
mt$classicFisher = gsub("<", "", mt$classicFisher) #some p-values have the "less than" sign ("<"), which causes the numeric column to be interpreted as character. 
mt$bonferroni = p.adjust(mt$classicFisher, method="bonferroni")
mt$bh = p.adjust(mt$classicFisher, method="BH")

if (!dir.exists(out_prefix)) {dir.create(out_prefix)}
out_file = paste0(c(ontology, "GO.results.tsv"), collapse = '.')
out_file = paste0(out_prefix, "/", out_file)

write.table(mt, file = out_file, sep='\t', quote = FALSE, row.names = FALSE)
}

calculate_go_enrichment('BP')
calculate_go_enrichment('MF')
calculate_go_enrichment('CC')





