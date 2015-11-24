##Download curated pathways from MSigDB: 
##http://software.broadinstitute.org/gsea/msigdb/download_file.jsp?filePath=/resources/msigdb/5.0/c2.cp.v5.0.symbols.gmt

library(data.table)
msigdb = "~/Desktop/RCAS/c2.cp.v5.0.symbols.gmt"

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

get_background = function(gene_lists){
  all_genes = c()
  for (l in gene_lists){
    all_genes = c(all_genes, l)
  }
  return(unique(all_genes))
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
    pval  = fisher.test(comparison, alternative = "two.sided")$p.value
    t_counts = c(t_counts, t)
    b_counts = c(b_counts, b)
    exp_vals = c(exp_vals, exp)
    pval_calc = c(pval_calc, pval)
  }
  results = cbind.data.frame(mynames, t_counts, rep(t_size, length(mynames)), b_counts, rep(b_size, length(mynames)), exp_vals, pval_calc)
  colnames(results) = c("list_name", "treatment_count", "treatment_size", "background_count", "background_size", "expected_in_treatment", "pval")
  results$bonferroni = p.adjust(results$pval, method = "bonferroni")
  results$BH = p.adjust(results$pval, method = "BH")
  return(data.table(results))
}


#parse lists of genes from MSigDB
gene_lists = parse_msigdb(msigdb)

#TODO substitute this with the list of all genes 
#get the background list of genes  
background = get_background(gene_lists)

#TODO substitute this with the list obtained from an experiment etc. (e.g. genes that overlap par-clip peaks)
#get the gene list of interest
treatment = gene_lists$KEGG_CITRATE_CYCLE_TCA_CYCLE

#count the number of times each pathway is associated to the gene list of interest
#compare the frequency of associations of different pathways to the treatment set versus the background set. 
#find pathways that are enriched in the treatment set. 
results = count_associations(treatment, background, gene_lists)
results = results[order(pval)]

#Write pathways with FDR < 0.001 to file
write.table(results[bonferroni < 0.001], file = "rcas.msigdb.results.tsv", quote = FALSE, row.names = FALSE, sep="\t")




