library(RCAS)
context("Functions to do gene set enrichment analysis")

data(gff)
data(queryRegions)
data(geneSets)
overlaps <- queryGff(queryRegions = queryRegions, gffData = gff)

results <- runGSEA(geneSetList = geneSets,
                     species = 'human',
                     backgroundGenes = unique(gff$gene_id),
                     targetedGenes = unique(overlaps$gene_id))

test_that("Gene set enrichment analysis using random gene sets", {
  expect_is(results, 'data.frame')
  expect_equal(colnames(results)[5], 'BH')
  expect_equal(colnames(results)[6], 'bonferroni')
  expect_equal(colnames(results)[7], 'foldEnrichment')
})

test_that("Creating biomart connections for different genome versions",{
  expect_is(getBioMartConnection('hsapiens'), 'Mart')
  expect_is(getBioMartConnection('mmusculus'), 'Mart')
  expect_is(getBioMartConnection('dmelanogaster'), 'Mart')
  expect_is(getBioMartConnection('celegans'), 'Mart')
  expect_error(getBioMartConnection('foo'))
  expect_error(getBioMartConnection('bar'))
})

genes <- c('2645','5232', '5230','5162','5160')
orth <- retrieveOrthologs(mart1 = getBioMartConnection('hsapiens'),
                          mart2 = getBioMartConnection('mmusculus'),
                          geneSet = genes)
test_that("Retrieving orthologs of a given set of genes", {
  expect_equal(length(orth$NCBI.gene.ID.1), 5)
})

