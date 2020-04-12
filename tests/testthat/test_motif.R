library(RCAS)
context("Functions to do motif enrichment analysis")

data(queryRegions)
peaks <- queryRegions[1:1000]

controls <- createControlRegions(queryRegions = peaks)
test_that("Generating control regions to compare with query regions", {
  expect_is(controls, 'GRanges')
  expect_equal(length(controls), length(peaks))
  expect_equal(sum(width(controls) == width(peaks)), length(peaks))
})

test_that("Extracting sequences of a GRanges object from a BSGenome object", {
  expect_is(extractSequences(peaks[1:5], 'hg19'), 'DNAStringSet')
  expect_equal(width(extractSequences(peaks[1:5], 'hg19')), width(peaks[1:5]))
  expect_error(extractSequences(peaks[1:5], 'foo'))
})

motifs <- runMotifDiscovery(queryRegions = queryRegions, genomeVersion = 'hg19', 
                            motifWidth = 4, maxMismatch = 0, 
                            motifN = 1, nCores = 1)
test_that("Motif analysis ... ",{
  expect_is(motifs, 'list')
})

test_that("Motif summary table", {
  expect_equal(length(colnames(getMotifSummaryTable(motifs))), 9)
  expect_equal(colnames(getMotifSummaryTable(motifs))[1], 'patterns')
})
