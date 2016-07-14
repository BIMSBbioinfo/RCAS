library(RCAS)
context("tests for main runReport function")

testBedFile <- system.file('extdata', 'testfile.bed', package = 'RCAS')

test_that("Testing runReport function", {
  expect_error(runReport(genomeVersion = 'foo'),
               regexp = "not a supported genome version.",
               all = FALSE)
  expect_error(runReport(genomeVersion = 'mm9'),
               regexp = "Test-data only works for human",
               all = FALSE)
  expect_error(runReport(queryFilePath = testBedFile, genomeVersion = 'mm9'),
               regexp = 'Test-data only works for human',
               all = FALSE)
})

