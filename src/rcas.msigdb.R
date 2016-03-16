##Download curated pathways from MSigDB: 
##http://software.broadinstitute.org/gsea/msigdb/download_file.jsp?filePath=/resources/msigdb/5.0/c2.cp.v5.0.symbols.gmt

suppressWarnings(suppressMessages(library('org.Hs.eg.db'))) #mouse: org.Mm.eg.db
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library('org.Ce.eg.db'))) 
suppressWarnings(suppressMessages(library('org.Mm.eg.db'))) 
suppressWarnings(suppressMessages(library('org.Dm.eg.db'))) 
suppressWarnings(suppressMessages(library('rtracklayer')))
  
#1. Collect arguments
args <- commandArgs(TRUE)

## Default setting when no arguments passed
if(length(args) < 4) {
  args <- c("--help")
}

help_command = "
rcas.msigdb.R: Gene Set Enrichment Analysis module of RCAS based on Gene Sets from MSigDB

Arguments:
--gmt=<path to MSigDB gmt file>  - e.g /home/buyar/projecst/RCAS/c2.cp.v5.0.symbols.gmt
--background_list=<background list of genes >  - e.g {sample}.background_genes.txt
--targeted_list=<peaks file in BED format>  - e.g {sample}.targeted_genes.txt
--out_prefix=<output prefix>   -e.g Hafner2010.hg19
--species=<species name>    -Choose between (human, fly, worm, mouse)
--help              - print this text

Example:
Rscript rcas.msigdb.R --gmt=c2.cp.v5.0.entrez.gmt --background_list={sample}.background_genes.txt--targeted_list={sample}.targeted_genes.txt --out_prefix=Hafner2010.hg19 --species=human"

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

if(!("gmt" %in% argsDF$V1)) {
  cat(help_command, "\n")
  stop("provide the path to gmt file")
}

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

msigdb = argsL$gmt
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
#2. Get the list of all protein-coding genes from the GTF file 

all_gene_ids = na.omit(scan(background_list, what = character()))
head(all_gene_ids)

cat("Read ",length(all_gene_ids),"genes from ",background_list,"\n")

# Get the query regions info
targeted_gene_ids = na.omit(scan(targeted_list, what = character()))
cat("Found ",length(targeted_gene_ids),"genes targeted by query regions\n")

#4. Look for enriched MSigDB gene sets for the genes found in the anot.bed compared to the background genes (found in gtf)
#source("http://www.bioconductor.org/biocLite.R")
#biocLite(c(org.Hs.eg.db"))

parse_msigdb = function(msigdb){
  data = readLines(msigdb)
  
  gene_lists = list()
  mynames = c()
  for (i in 1:length(data)){
    info = unlist(strsplit(data[i], "\t"))
    gene_lists = c(gene_lists, list(paste(info[3:length(info)], sep= "\t")))
    mynames = c(mynames, info[1])
  }
  names(gene_lists) = mynames	
  
  return(gene_lists)
}

count_associations = function(treatment, background, gene_lists){
  i = 0
  mynames = names(gene_lists)
  t_size = length(treatment)
  b_size = length(background)
  t_counts = c()
  b_counts = c()
  exp_vals = c()
  pval_calc = c()
  for (l in gene_lists){
    i = i + 1
    name = mynames[i]
    t = sum(treatment %in% l)
    b = sum(background %in% l)
    exp = t_size * (b/b_size)
    comparison =  matrix(c(t, b, t_size - t, b_size - b), nrow = 2, dimnames = list(c("treatment", "background"), c("found", "not_found")))
    pval  = gsub(pattern = '<', replacement = '', x = fisher.test(comparison, alternative = "greater")$p.value)
    
    t_counts = c(t_counts, t)
    b_counts = c(b_counts, b)
    exp_vals = c(exp_vals, exp)
    pval_calc = c(pval_calc, pval)
  }
  #for debugging#
  #results = cbind.data.frame(mynames, t_counts, rep(t_size, length(mynames)), b_counts, rep(b_size, length(mynames)), exp_vals, pval_calc)
  #colnames(results) = c("list_name", "treatment_count", "treatment_size", "background_count", "background_size", "expected_in_treatment", "pval")
  
  results = cbind.data.frame(mynames, t_counts, rep(t_size, length(mynames)), round(exp_vals, 1), pval_calc)
  colnames(results) = c("list_name", "treatment_count", "treatment_size", "expected_in_treatment", "pval")
  
  results$bonferroni = format.pval(p.adjust(results$pval, method = "bonferroni"))
  results$BH = format.pval(p.adjust(results$pval, method = "BH"))
  results$pval = format.pval(results$pval)
  return(data.table(results))
}

get_unique_items_from_list = function(mylist){
  
  items = c()
  for (l in mylist)
  {
    items = c(items, unique(l))
  }
  return(items)
}

#parse lists of genes from MSigDB 
gene_lists = parse_msigdb(msigdb)
cat("got the gene lists from",msigdb,"\n")

#map ENSEMBL gene ids to ENtrez Gene Ids
if (species == 'human'){
  mapping_dict = org.Hs.egENSEMBL2EG
}else if (species == 'mouse'){
  mapping_dict = org.Mm.egENSEMBL2EG
}else if (species == 'fly'){
  mapping_dict = org.Dm.egENSEMBL2EG
}else if (species == 'worm'){
  mapping_dict = org.Ce.egENSEMBL2EG
}
ens2eg <- as.list(mapping_dict) #mapping dictionary from ENSEMBL to ENTREZ 
treatment = get_unique_items_from_list(ens2eg[gsub("\\..*$", "", targeted_gene_ids)])
background = get_unique_items_from_list(ens2eg[gsub("\\..*$", "", all_gene_ids)])

cat("mapped ",length(targeted_gene_ids),"from ensembl to",length(treatment),"gene ids from entrez for treatment", "\n")
cat("mapped ",length(all_gene_ids)," from ensembl to",length(background),"gene ids from entrez for background", "\n")

#count the number of times each pathway is associated to the gene list of interest
#compare the frequency of associations of different pathways to the treatment set versus the background set. 
#find pathways that are enriched in the treatment set. 
results = count_associations(treatment = treatment, background = background, gene_lists)
results = results[order(pval)]

#Write pathways with FDR < 0.001 to file
#write.table(results[bonferroni < 0.001], file = out_file, quote = FALSE, row.names = FALSE, sep="\t")
out_file = paste0(out_prefix,'.msigdb_results.tsv')
write.table(results, file = out_file, quote = FALSE, row.names = FALSE, sep="\t")




