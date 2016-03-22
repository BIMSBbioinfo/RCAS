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
