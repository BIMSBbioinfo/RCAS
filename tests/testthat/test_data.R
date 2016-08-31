library(RCAS)
context("built-in data: gff, queryRegions, msigDB")

data(gff)
test_that("gff is a GRanges object containing test data of dimensions", {
  expect_is(object = gff, 'GRanges')
  expect_equal(length(gff), 238010)
  expect_equal(ncol(GenomicRanges::mcols(gff)), 16)
})

data(queryRegions)
test_that("queryRegions is a GRanges object from a HITS-CLIP dataset", {
  expect_is(object = queryRegions, 'GRanges')
  expect_equal(length(queryRegions), 10000)
  expect_equal(round(mean(GenomicRanges::score(queryRegions))), 45)
})

data(geneSets)
test_that("geneSets is a list of vectors containing gene sets", {
  expect_is(object = geneSets, 'list')
  expect_equal(length(geneSets), 100)
  expect_match(names(msigDB)[1], 'randomGeneSet1')
})
