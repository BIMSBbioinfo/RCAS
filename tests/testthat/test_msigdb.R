library(RCAS)
context("Functions to do MSIGDB gene set enrichment analysis")

data(gff)
data(queryRegions)
data(msigDB)
overlaps <- queryGff(queryRegions = queryRegions, gff = gff)

results <- runMSIGDB(msigDB = msigDB,
                     species = 'human',
                     backgroundGenes = unique(gff$gene_id),
                     targetedGenes = unique(overlaps$gene_id))

test_that("Gene set enrichment analysis using MSIGDB gene sets", {
  expect_is(results, 'data.frame')
  expect_equal(colnames(results)[5], 'BH')
  expect_equal(colnames(results)[6], 'bonferroni')
})

test_that("Creating biomart connections for different genome versions",{
  expect_is(getBioMartConnection(genomeVersion = 'hg19'), 'Mart')
  expect_is(getBioMartConnection(genomeVersion = 'dm3'), 'Mart')
  expect_is(getBioMartConnection(genomeVersion = 'ce10'), 'Mart')
  expect_is(getBioMartConnection(genomeVersion = 'mm9'), 'Mart')
  expect_error(getBioMartConnection(genomeVersion = 'mm10'))
  expect_error(getBioMartConnection(genomeVersion = 'bar'))
  expect_error(getBioMartConnection(genomeVersion = 'foo'))
})

genes <- c('2645','5232', '5230','5162','5160')
orth <- retrieveOrthologs(mart1 = getBioMartConnection('hg19'),
                          mart2 = getBioMartConnection('mm9'),
                          geneSet = genes)
test_that("Retrieving orthologs of a given set of genes", {
  expect_equal(length(orth$EntrezGene.ID.1), 5)
})

