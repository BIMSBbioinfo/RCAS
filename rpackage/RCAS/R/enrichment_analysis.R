#functions for GO term and MSIGDB gene set enrichment analyses

runTopGO <- function (ontology = 'BP', species = 'human', backgroundGenes, targetedGenes) {

  category = data.frame(geneId = backgroundGenes)
  category$targeted = 0
  category[category$geneId %in% targetedGenes,]$targeted = 1
  category$geneId = gsub("\\..*$", "", category$geneId)
  allGenes = category$targeted
  names(allGenes) = category$geneId

  if (species == 'human') {
    speciesOrgdb = 'org.Hs.eg.db'
  } else if (species == 'mouse') {
    speciesOrgdb = 'org.Mm.eg.db'
  } else if (species == 'fly') {
    speciesOrgdb = 'org.Dm.eg.db'
  } else if (species == 'worm') {
    speciesOrgdb = 'org.Ce.eg.db'
  } else {
    stop ("Cannot do GO term analysis for species except human, worm, fly, and mouse\n")
  }

  GOdata <- new("topGOdata", ontology = ontology, allGenes = allGenes, geneSel = function (x){ return(x == 1) }, nodeSize=20, description = "Test", annot = annFUN.org, mapping=speciesOrgdb, ID="Ensembl")

  resultFisher <- topGO::runTest(GOdata, algorithm = "classic", statistic = "fisher")
  goResults = topGO::GenTable(GOdata, classicFisher = resultFisher, topNodes = length(usedGO(GOdata)))

  #Do multiple-testing correction
  goResults$classicFisher = gsub("<", "", goResults$classicFisher) #some p-values have the "less than" sign ("<"), which causes the numeric column to be interpreted as character.
  goResults$bonferroni = p.adjust(goResults$classicFisher, method="bonferroni")
  goResults$bh = p.adjust(goResults$classicFisher, method="BH")
  return(goResults)
}

parseMsigdb = function(filePath){
  data = readLines(filePath)
  geneLists = list()
  geneListNames = c()
  for (i in 1:length(data)){
    info = unlist(strsplit(data[i], "\t"))
    geneLists = c(geneLists, list(paste(info[3:length(info)], sep= "\t")))
    geneListNames = c(geneListNames, info[1])
  }
  names(geneLists) = geneListNames
  return(geneLists)
}

#given two biomart connections and a set of entrez gene identifiers; retrieve orthologs between mart1 and mart2 for the given list of genes
retrieveOrthologs = function(mart1, mart2, geneSet){
  biomaRt::getLDS(attributes = c("entrezgene"),
         filters = "entrezgene", values = geneSet, mart = mart1,
         attributesL = c("entrezgene"), martL = mart2)
}

getBioMartConnection <- function (genomeVersion) {
  if (genomeVersion == 'hg19') {
    mart <- biomaRt::useMart(biomart='ENSEMBL_MART_ENSEMBL', host='feb2014.archive.ensembl.org', dataset = "hsapiens_gene_ensembl")
  } else if (genomeVersion == 'mm9') {
    mart <- biomaRt::useMart(biomart='ENSEMBL_MART_ENSEMBL', host='may2012.archive.ensembl.org', dataset = "mmusculus_gene_ensembl")
  } else if (genomeVersion == 'dm3') {
    mart <- biomaRt::useMart(biomart='ENSEMBL_MART_ENSEMBL', host='dec2014.archive.ensembl.org', dataset = "dmelanogaster_gene_ensembl")
  } else if (genomeVersion == 'ce6') {
    mart <- biomaRt::useMart(biomart='ENSEMBL_MART_ENSEMBL', host='may2009.archive.ensembl.org', dataset = "celegans_gene_ensembl")
  } else if (genomeVersion == 'ce10') {
    mart <- biomaRt::useMart(biomart='ENSEMBL_MART_ENSEMBL', host='jan2013.archive.ensembl.org', dataset = "celegans_gene_ensembl")
  } else {
    stop ("Cannot create a BioMart connection for genome versions except hg19, ce6, ce10, dm3, and mm9\n")
  }
  return (mart)
}

#' @export
createOrthologousMsigdbDataset <- function(refMsigdbFilePath, refGenomeVersion, targetGenomeVersion){
  #parse lists of genes from MSigDB
  referenceGeneSets = parseMsigdb(refMsigdbFilePath)
  cat("got",length(names(referenceGeneSets)),"gene lists from",refMsigdbFilePath,"\n")

  #define all available reference genes in all the given gene sets
  refGenes = unique(paste(unlist(referenceGeneSets)))
  cat('Retrieved',length(refGenes),'genes for reference genome in the MSIGDB input file\n')

  # to be able to retrieve orthologs of the genes in the input gene sets - use biomaRt to get ortholog data from Ensembl
  refMart <- getBioMartConnection (genomeVersion = refGenomeVersion)
  targetMart <- getBioMartConnection (genomeVersion = targetGenomeVersion)
  cat('Created the database connections to biomart for',refGenomeVersion,'and',targetGenomeVersion,'\n')

  #retrieve orthologs lists for all human genes found in the MSIGDB gene sets
  orthologs = retrieveOrthologs(mart1 = refMart, mart2 = targetMart, refGenes)
  cat('Retrieved',nrow(orthologs),'orthologous relations between',length(unique(orthologs[,1])),"genes from",refGenomeVersion,'and',length(unique(orthologs[,2])),"genes from",targetGenomeVersion,"\n")

  orthMsigdbDataset = list()
  counts1 = c()
  counts2 = c()
  for (i in 1:length(referenceGeneSets)){
    refSet = unique(as.numeric(paste(unlist(referenceGeneSets[i]))))
    orthologSet = unique(orthologs[orthologs$EntrezGene.ID %in% refSet,]$EntrezGene.ID.1)
    refSetName = names(referenceGeneSets)[i]
    orthMsigdbDataset[[refSetName]] = orthologSet
    counts1 = c(counts1, length(refSet))
    counts2 = c(counts2, length(orthologSet))
  }
  return(orthMsigdbDataset)
}

runMSIGDB <- function (msigDB, species = 'human', backgroundGenes, targetedGenes) {
  #map ENSEMBL gene ids to Entrez Gene Ids
  if (species == 'human'){
    mappingDict = org.Hs.egENSEMBL2EG
  }else if (species == 'mouse'){
    mappingDict = org.Mm.egENSEMBL2EG
  }else if (species == 'fly'){
    mappingDict = org.Dm.egENSEMBL2EG
  }else if (species == 'worm'){
    mappingDict = org.Ce.egENSEMBL2EG
  }
  ens2eg <- as.list(mappingDict) #mapping dictionary from ENSEMBL to ENTREZ
  treatment = get_unique_items_from_list(ens2eg[gsub("\\..*$", "", targeted_gene_ids)])
  background = get_unique_items_from_list(ens2eg[gsub("\\..*$", "", all_gene_ids)])
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
    pval  = fisher.test(comparison, alternative = "greater")$p.value

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




