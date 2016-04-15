library(RCAS)
context("tests for main runReport function")

test_that("Testing runReport function", {
  expect_error(runReport(genomeVersion = 'foo'), "foo is not a supported genome version.")
  expect_error(runReport(genomeVersion = 'mm9'), "\n  Test-data only works for human.\n           Please provide a queryFilePath to input in BED format \n\n")
  expect_error(runReport(queryFilePath = 'foo', genomeVersion = 'mm9'), '\n  Test-data only works for human.\n           Please provide a gffFilePath to input in GTF format \n\n')
  expect_error(runReport(queryFilePath = 'foo', gffFilePath = 'bar', genomeVersion='mm9'), ' \n  Test-data only works for human.\n           Please provide a gene set dataset with ENTREZ gene ids\n            downloaded from MSIGDB database or\n            set msigdbAnalysis option to FALSE \n\n')
})

