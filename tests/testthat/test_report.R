library(RCAS)
context("tests for main runReport function")

test_that("Testing runReport function", {
  expect_error(runReport(genomeVersion = 'foo'),
               regexp = "Can't find a genome sequence")
  expect_error(runReport(genomeVersion = '9'),
               regexp = "match the genome version to an unambigous genome sequence")
})

