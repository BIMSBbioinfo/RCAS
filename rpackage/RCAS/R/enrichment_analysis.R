#functions for GO term and MSIGDB gene set enrichment analyses

runTopGO <- function (ontology = 'BP', species = 'human', backgroundGenes, targetedGenes) {

  category <- data.frame(geneId = backgroundGenes)
  category$targeted <- 0
  category[category$geneId %in% targetedGenes,]$targeted <- 1
  category$geneId <- gsub("\\..*$", "", category$geneId)
  allGenes <- category$targeted
  names(allGenes) <- category$geneId

  if (species == 'human') {
    speciesOrgdb <- 'org.Hs.eg.db'
  } else if (species == 'mouse') {
    speciesOrgdb <- 'org.Mm.eg.db'
  } else if (species == 'fly') {
    speciesOrgdb <- 'org.Dm.eg.db'
  } else if (species == 'worm') {
    speciesOrgdb <- 'org.Ce.eg.db'
  } else {
    stop ("Cannot do GO term analysis for species except human, worm, fly, and mouse\n")
  }

  GOdata <- new("topGOdata", ontology = ontology, allGenes = allGenes, geneSel = function (x){ return(x == 1) }, nodeSize=20, description = "Test", annot = annFUN.org, mapping=speciesOrgdb, ID="Ensembl")

  resultFisher <- topGO::runTest(GOdata, algorithm = "classic", statistic = "fisher")
  goResults <- topGO::GenTable(GOdata, classicFisher = resultFisher, topNodes = length(usedGO(GOdata)))

  #Do multiple-testing correction
  goResults$classicFisher <- gsub("<", "", goResults$classicFisher) #some p-values have the "less than" sign ("<"), which causes the numeric column to be interpreted as character.
  goResults$bonferroni <- p.adjust(goResults$classicFisher, method="bonferroni")
  goResults$bh <- p.adjust(goResults$classicFisher, method="BH")
  return(goResults)
}

#' @export
parseMsigdb <- function(filePath){
  data <- readLines(filePath)
  geneLists <- list()
  geneListNames <- c()
  for (i in 1:length(data)){
    info <- unlist(strsplit(data[i], "\t"))
    geneLists <- c(geneLists, list(paste(info[3:length(info)], sep= "\t")))
    geneListNames <- c(geneListNames, info[1])
  }
  names(geneLists) <- geneListNames
  return(geneLists)
}

#given two biomart connections and a set of entrez gene identifiers; retrieve orthologs between mart1 and mart2 for the given list of genes
retrieveOrthologs <- function(mart1, mart2, geneSet){
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
  referenceGeneSets <- parseMsigdb(refMsigdbFilePath)
  cat("got",length(names(referenceGeneSets)),"gene lists from",refMsigdbFilePath,"\n")

  #define all available reference genes in all the given gene sets
  refGenes <- unique(paste(unlist(referenceGeneSets)))
  cat('Retrieved',length(refGenes),'genes for reference genome in the MSIGDB input file\n')

  # to be able to retrieve orthologs of the genes in the input gene sets - use biomaRt to get ortholog data from Ensembl
  refMart <- getBioMartConnection (genomeVersion = refGenomeVersion)
  targetMart <- getBioMartConnection (genomeVersion = targetGenomeVersion)
  cat('Created the database connections to biomart for',refGenomeVersion,'and',targetGenomeVersion,'\n')

  #retrieve orthologs lists for all human genes found in the MSIGDB gene sets
  orthologs <- retrieveOrthologs(mart1 = refMart, mart2 = targetMart, refGenes)
  cat('Retrieved',nrow(orthologs),'orthologous relations between',length(unique(orthologs[,1])),"genes from",refGenomeVersion,'and',length(unique(orthologs[,2])),"genes from",targetGenomeVersion,"\n")

  orthMsigdbDataset <- list()
  for (i in 1:length(referenceGeneSets)){
    refSet <- unique(as.numeric(paste(unlist(referenceGeneSets[i]))))
    orthologSet <- unique(orthologs[orthologs$EntrezGene.ID %in% refSet,]$EntrezGene.ID.1)
    refSetName <- names(referenceGeneSets)[i]
    orthMsigdbDataset[[refSetName]] <- orthologSet
  }
  return(orthMsigdbDataset)
}

calculateEnrichment <- function (targetedGenes, backgroundGenes, geneSet) {

  #we need to know the number of all genes in the compared sets
  tSize <- length(targetedGenes)
  bSize <- length(backgroundGenes)

  #we need to know how many genes from targeted and background lists are members of the given gene set
  t <- sum(targetedGenes %in% geneSet)
  b <- sum(backgroundGenes %in% geneSet)

  #how many genes do we expect find in the given gene set, based on the number of genes in targeted gene sets relative to the background
  exp <- round(tSize * (b/bSize), 1)

  contingencyTable <- matrix(c(t, b, tSize - t, bSize - b), nrow = 2, dimnames = list(c("treatment", "background"), c("found", "not_found")))
  pval  <- fisher.test(contingencyTable, alternative = "greater")$p.value

  result <- data.frame('treatment' = t, 'treatmentSize' = tSize, 'expectedInTreatment' = exp, 'fisherPVal' = pval)
  return (result)
}

#' @export
runMSIGDB <- function (msigDB, species = 'human', backgroundGenes, targetedGenes) {
  #map ENSEMBL gene ids to Entrez Gene Ids
  if (species == 'human'){
    mappingDict <- org.Hs.egENSEMBL2EG
  }else if (species == 'mouse'){
    mappingDict <- org.Mm.egENSEMBL2EG
  }else if (species == 'fly'){
    mappingDict <- org.Dm.egENSEMBL2EG
  }else if (species == 'worm'){
    mappingDict <- org.Ce.egENSEMBL2EG
  }
  ens2eg <- AnnotationDbi::as.list(mappingDict) #mapping dictionary from ENSEMBL to ENTREZ
  #map gene ids from ENSEMBL to ENTREZ
  targetedGenes <- unique(paste(unlist(ens2eg[targetedGenes])))
  backgroundGenes <- unique(paste(unlist(ens2eg[backgroundGenes])))

  results <- lapply(X = msigDB,
                   FUN = function(x) {
                     calculateEnrichment(targetedGenes=targetedGenes,
                                         backgroundGenes=backgroundGenes,
                                         geneSet = x )})
  results <- do.call("rbind", results)
  results$BH <- p.adjust(results$fisherPVal, method = "BH")
  results$bonferroni <- p.adjust(results$fisherPVal, method = "bonferroni")
  return(results[order(results$fisherPVal),])
}




