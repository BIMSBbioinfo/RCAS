library(RCAS)
suppressMessages(library(topGO))
context("Functions to do GO term enrichment analysis")

data(gff)
data(queryRegions)
test_that("Importing test data", {
  expect_is(gff, 'GRanges')
  expect_is(queryRegions, 'GRanges')
})

overlaps <- queryGff(queryRegions = queryRegions, gffData = gff)
test_that("Testing queryGff function", {
  expect_is(overlaps, 'GRanges')
  expect_equal(length(overlaps), 12200)
  expect_equal(length(colnames(GenomicRanges::mcols(overlaps))), 17)
})


goResults <- runTopGO(ontology = 'CC',
                      species = 'human',
                      backgroundGenes = unique(gff$gene_id),
                      targetedGenes = unique(overlaps$gene_id))
test_that("GO term analysis for cellular compartments", {
  expect_is(goResults, 'data.frame')
  expect_equal(colnames(goResults)[7], 'bonferroni')
  expect_equal(colnames(goResults)[8], 'bh')
})

