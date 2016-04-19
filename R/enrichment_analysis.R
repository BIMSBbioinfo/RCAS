.datatable.aware=TRUE

#' runTopGO
#'
#' A wrapper function to facilitate GO term enrichment analysis using
#' \code{topGO} package
#'
#' @param ontology A character string denoting which type of GO ontology to use.
#'   Options are BP (biological processes), MF (molecular functions) and CC
#'   (cellular compartments).
#' @param species A character string denoting which species is under analysis.
#'   Options are 'human', 'mouse', 'fly' and 'worm'.
#' @param backgroundGenes A vector of Ensembl gene ids that serve as background
#'   set of genes for GO term enrichment. In the context of RCAS, this should be
#'   the whole set of genes found in an input GTF data.
#' @param targetedGenes A vector of Ensembl gene ids that serve as the set for
#'   which GO term enrichment should be carried out. In the context of RCAS,
#'   this should be the set of genes that overlap with the query regions in an
#'   input BED file.
#'
#' @return A data.frame object containing enriched GO terms and associated
#'   statistics
#'
#' @examples
#' \dontrun{
#' #get all genes from the gff data
#' backgroundGenes <- unique(gff$gene_id)
#' #get genes that overlap query regions
#' overlaps <- queryGff(queryRegions, gff)
#' targetedGenes <- unique(overlaps$gene_id)
#' #run TopGO
#' goResults = runTopGO( ontology = 'BP',
#'                       species = 'human',
#'                       backgroundGenes = backgroundGenes,
#'                       targetedGenes = targetedGenes)
#'                       }
#' @import org.Hs.eg.db
#' @importFrom topGO runTest
#' @importFrom topGO GenTable
#' @importFrom topGO usedGO
#' @importFrom topGO annFUN.org
#' @importFrom stats p.adjust
#' @importFrom stats fisher.test
#' @export
runTopGO <- function (ontology = 'BP',
                      species = 'human',
                      backgroundGenes,
                      targetedGenes) {

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
    stop ("Cannot do GO term analysis for species except:
          human, worm, fly, and mouse\n")
  }

  GOdata <- new("topGOdata",
                ontology = ontology,
                allGenes = allGenes,
                geneSelectionFun = function (x){ return(x == 1) },
                nodeSize = 20,
                description = "Test",
                annot = annFUN.org,
                mapping = speciesOrgdb,
                ID = "Ensembl")

  resultFisher <- topGO::runTest(object = GOdata,
                                 algorithm = "classic",
                                 statistic = "fisher")
  goResults <- topGO::GenTable(object = GOdata,
                               classicFisher = resultFisher,
                               topNodes = length(topGO::usedGO(GOdata)))

  #Do multiple-testing correction
  #some p-values have the "less than" sign ("<"),
  #which causes the numeric column to be interpreted as character.
  goResults$classicFisher <- gsub("<", "", goResults$classicFisher)
  goResults$bonferroni <- stats::p.adjust(p = goResults$classicFisher,
                                          method = "bonferroni")
  goResults$bh <- stats::p.adjust(goResults$classicFisher, method="BH")
  return(goResults)
}

#' parseMsigdb
#'
#' A function to import gene sets downloaded from the Molecular Signatures
#' Database (MSIGDB)
#'
#' @param filePath Path to a file containing gene sets from MSIGDB. The gene ids
#'   must be in Entrez format.
#'
#' @return A list of vectors where each vector consists of a set of Entrez gene
#'   ids
#'
#' @examples
#' #First Download gene sets (with Entrez Ids) from MSIGDB database
#' #from \url{http://software.broadinstitute.org/gsea/msigdb/collections.jsp#C2}
#'
#' input <- system.file('msigdb_test.gmt', package='RCAS')
#' msigDB <- parseMsigdb (filePath = input)
#'
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

#' Print MSIGDB Dataset to a file
#' This function is used to print a MSIGDB dataset into a file. Mostly
#' useful when human data is mapped to another species, and that mapping
#' is required to run the report.
#' @param dataset A list of vectors containing gene sets from MSIGDB
#' @param outputFilename A character string that denotes the output file name
#' @return A text file printed to the current directory
#' @examples
#' data(msigDB)
#' printMsigdbDataset(msigDB, 'output.gmt')
#'
#' @export
printMsigdbDataset = function(dataset, outputFilename){
  if (file.exists(outputFilename)){
    stop('A filename with the same name exists.
        Delete the file first or choose a different file')
  }
  sink(file = outputFilename, append = TRUE)
  for (i in 1:length(dataset)){
    geneSetName = names(dataset)[i]
    geneSet = paste(unique(paste(unlist(dataset[i]))), collapse = '\t')
    cat(geneSetName, 'no_http', geneSet, sep='\t', "\n")
  }
  sink()
}

#' retrieveOrthologs
#'
#' Given two biomart connections and a set of entrez gene identifiers; retrieve
#' orthologs between mart1 and mart2 for the given list of genes
#'
#' @param mart1 An Ensembl biomart connection for reference species created
#'   using the \code{biomaRt::useMart()} function
#' @param mart2 An Ensembl biomart connection for target species created using
#'   the \code{biomaRt::useMart()} function
#' @param geneSet A vector of Entrez gene ids from a reference species (should
#'   be available at the biomart object, mart1)
#' @return A data.frame object containing a mapping of orthologouse genes from
#'   two mart objects
#'
#' @examples
#' mart1_hg19 <- biomaRt::useMart(biomart = 'ENSEMBL_MART_ENSEMBL',
#'                                   host = 'feb2014.archive.ensembl.org',
#'                                dataset = "hsapiens_gene_ensembl")
#' mart2_mm9 <- biomaRt::useMart(biomart = 'ENSEMBL_MART_ENSEMBL',
#'                                  host = 'may2012.archive.ensembl.org',
#'                               dataset = "mmusculus_gene_ensembl")
#' genes <- c('2645','5232', '5230','5162','5160')
#' orthologs <- retrieveOrthologs( mart1 = mart1_hg19,
#'                                 mart2 = mart2_mm9,
#'                                 geneSet = genes)
#' @importFrom biomaRt getLDS
#' @export
retrieveOrthologs <- function(mart1, mart2, geneSet){
  biomaRt::getLDS( attributes = c("entrezgene"),
                   filters = "entrezgene",
                   values = geneSet,
                   mart = mart1,
                   attributesL = c("entrezgene"),
                   martL = mart2)
}

#' @importFrom biomaRt useMart
getBioMartConnection <- function (genomeVersion) {
  if (genomeVersion == 'hg19') {
    mart <- biomaRt::useMart( biomart = 'ENSEMBL_MART_ENSEMBL',
                              host = 'feb2014.archive.ensembl.org',
                              dataset = "hsapiens_gene_ensembl")
  } else if (genomeVersion == 'mm9') {
    mart <- biomaRt::useMart( biomart = 'ENSEMBL_MART_ENSEMBL',
                              host = 'may2012.archive.ensembl.org',
                              dataset = "mmusculus_gene_ensembl")
  } else if (genomeVersion == 'dm3') {
    mart <- biomaRt::useMart( biomart = 'ENSEMBL_MART_ENSEMBL',
                              host = 'dec2014.archive.ensembl.org',
                              dataset = "dmelanogaster_gene_ensembl")
  } else if (genomeVersion == 'ce10') {
    mart <- biomaRt::useMart( biomart = 'ENSEMBL_MART_ENSEMBL',
                              host = 'jan2013.archive.ensembl.org',
                              dataset = "celegans_gene_ensembl")
  } else {
    stop ("Cannot create a BioMart connection for genome versions except:
          hg19, ce10, dm3, and mm9\n")
  }
  return (mart)
}

#' createOrthologousMsigdbDataset
#'
#' MSIGDB curates gene sets for human but not for other species. This function
#' is used to create such gene sets for other species such as mouse, fly, and
#' worm via orthologous relationships to human genes.
#'
#' @param refMsigdbFilePath Path to a file containing gene sets from MSIGDB. The
#'   gene ids must be in Entrez format.
#' @param refGenomeVersion Genome version of a reference species. (default:hg19)
#' @param targetGenomeVersion Genome version of a target species. Available
#'   options are mm9, dm3, and ce10
#'
#' @return A list of vectors where each vector consists of a set of Entrez gene
#'   ids
#'
#' @examples
#' #First Download gene sets (with Entrez Ids) from MSIGDB database
#' #from \url{http://software.broadinstitute.org/gsea/msigdb/collections.jsp#C2}
#'
#' #Map the gene sets to a target genome (supported genomes: mm9, dm3, or ce10)
#' \dontrun{
#' createOrthologousMsigdbDataset(refMsigdbFilePath = 'msigdb.entrez.txt',
#'                                refGenomeVersion = 'hg19',
#'                                targetGenomeVersion = 'mm9')
#'                                }
#' @export
createOrthologousMsigdbDataset <- function(refMsigdbFilePath,
                                           refGenomeVersion = 'hg19',
                                           targetGenomeVersion) {
  #parse lists of genes from MSigDB
  referenceGeneSets <- parseMsigdb(filePath = refMsigdbFilePath)
  cat("got", length(names(referenceGeneSets)),
      "gene lists from", refMsigdbFilePath, "\n")

  #define all available reference genes in all the given gene sets
  refGenes <- unique(paste(unlist(referenceGeneSets)))
  cat('Retrieved', length(refGenes),
      'genes for reference genome in the MSIGDB input file\n')

  # to be able to retrieve orthologs of the genes in the input gene sets - use
  # biomaRt to get ortholog data from Ensembl
  refMart <- getBioMartConnection (genomeVersion = refGenomeVersion)
  targetMart <- getBioMartConnection (genomeVersion = targetGenomeVersion)
  cat('Created the database connections to biomart for',
      refGenomeVersion, 'and', targetGenomeVersion,'\n')

  #retrieve orthologs lists for all human genes found in the MSIGDB gene sets
  orthologs <- retrieveOrthologs(mart1 = refMart, mart2 = targetMart, refGenes)
  cat('Retrieved', nrow(orthologs),
      'orthologous relations between',
      length(unique(orthologs[,1])),
      "genes from", refGenomeVersion,
      'and', length(unique(orthologs[,2])),
      "genes from", targetGenomeVersion,"\n")

  orthMsigdbDataset <- list()
  for (i in 1:length(referenceGeneSets)){
    refSet <- unique(as.numeric(paste(unlist(referenceGeneSets[i]))))
    s <- orthologs$EntrezGene.ID %in% refSet
    orthologSet <- unique(orthologs[s,]$EntrezGene.ID.1)
    refSetName <- names(referenceGeneSets)[i]
    orthMsigdbDataset[[refSetName]] <- orthologSet
  }
  return(orthMsigdbDataset)
}

#' @importFrom stats fisher.test
calculateEnrichment <- function (targetedGenes, backgroundGenes, geneSet) {

  #we need to know the number of all genes in the compared sets
  tSize <- length(targetedGenes)
  bSize <- length(backgroundGenes)

  #we need to know how many genes from targeted and background lists are members
  #of the given gene set
  t <- sum(targetedGenes %in% geneSet)
  b <- sum(backgroundGenes %in% geneSet)

  #how many genes do we expect find in the given gene set, based on the number
  #of genes in targeted gene sets relative to the background
  exp <- round(tSize * (b/bSize), 1)

  contingencyTable <- matrix(c(t, b, tSize - t, bSize - b), nrow = 2,
                             dimnames = list(c("treatment", "background"),
                                             c("found", "not_found")))
  pval  <- stats::fisher.test(contingencyTable, alternative = "greater")$p.value

  result <- data.frame('treatment' = t,
                       'treatmentSize' = tSize,
                       'expectedInTreatment' = exp,
                       'fisherPVal' = pval)
  return (result)
}

#' runMSIGDB
#'
#' MSIGDB is a database of curated gene sets. This function is used to
#' facilitate gene set enrichment for genes that are found to overlap query
#' regions.
#'
#' @param msigDB A list of vectors where each vector consists of a set of Entrez
#'   gene ids returned by \code{parseMsigdb} function
#' @param species A character string denoting the species under analysis.
#'   Options are 'human', 'mouse', 'fly' and 'worm'.
#' @param backgroundGenes A vector of Ensembl gene ids that serve as background
#'   set of genes for GO term enrichment. In the context of RCAS, this should be
#'   the whole set of genes found in an input GTF file.
#' @param targetedGenes A vector of Ensembl gene ids that serve as the set for
#'   which GO term enrichment should be carried out. In the context of RCAS,
#'   this should be the set of genes that overlap with the query regions in an
#'   input BED file.
#' @return A data.frame object containing enriched MSIGDB gene sets and
#'   associated statistics
#'
#' @examples
#'
#' #load test data
#' data(msigDB)
#' data(gff)
#' data(queryRegions)
#' #get all genes from the gff data
#' backgroundGenes <- unique(gff$gene_id)
#' #get genes that overlap query regions
#' overlaps <- queryGff(queryRegions, gff)
#' targetedGenes <- unique(overlaps$gene_id)
#' runMSIGDB(msigDB = msigDB,
#'           species = 'human',
#'           backgroundGenes = backgroundGenes,
#'           targetedGenes = targetedGenes)
#' @importFrom org.Hs.eg.db org.Hs.egENSEMBL2EG
#' @importFrom AnnotationDbi as.list
#' @importFrom stats p.adjust
#' @export
runMSIGDB <- function (msigDB,
                       species = 'human',
                       backgroundGenes,
                       targetedGenes) {
  #map ENSEMBL gene ids to Entrez Gene Ids
  if (species == 'human'){
    mappingDict <- org.Hs.eg.db::org.Hs.egENSEMBL2EG
  }else if (species == 'mouse'){
    mappingDict <- org.Mm.eg.db::org.Mm.egENSEMBL2EG
  }else if (species == 'fly'){
    mappingDict <- org.Dm.eg.db::org.Dm.egENSEMBL2EG
  }else if (species == 'worm'){
    mappingDict <- org.Ce.eg.db::org.Ce.egENSEMBL2EG
  }
  #mapping dictionary from ENSEMBL to ENTREZ
  ens2eg <- AnnotationDbi::as.list(mappingDict)

  #map gene ids from ENSEMBL to ENTREZ
  targetedGenes <- unique(paste(unlist(ens2eg[targetedGenes])))
  backgroundGenes <- unique(paste(unlist(ens2eg[backgroundGenes])))

  results <- lapply(X = msigDB,
                   FUN = function(x) {
                     calculateEnrichment(targetedGenes=targetedGenes,
                                         backgroundGenes=backgroundGenes,
                                         geneSet = x )})
  results <- do.call("rbind", results)
  results$BH <- stats::p.adjust(results$fisherPVal, method = "BH")
  results$bonferroni <- stats::p.adjust(p = results$fisherPVal,
                                        method = "bonferroni")
  return(results[order(results$fisherPVal),])
}




