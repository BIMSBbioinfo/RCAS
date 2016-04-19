library(RCAS)
context("Functions to get txdbFeatures")

data(gff)
txdb <- GenomicFeatures::makeTxDbFromGRanges(gff)
test_that("Importing txdb", {
  expect_is(txdb, 'TxDb')
})

features <- getTxdbFeatures(txdb)
test_that("Parsing features from txdb", {
  expect_is(features, 'list')
  expect_equal(length(names(features)), 8)
  expect_match(names(features)[1], "transcripts")
  expect_match(names(features)[8], "threeUTRs")
})

features <- getTxdbFeaturesFromGff(gff = gff)
test_that("Parsing features from GFF", {
  expect_is(features, 'list')
  expect_equal(length(names(features)), 8)
  expect_match(names(features)[1], "transcripts")
  expect_match(names(features)[8], "threeUTRs")
})

data(queryRegions)
mytable <- getTargetedGenesTable(queryRegions = queryRegions,
                                 txdbFeatures = features)
test_that("Getting table of targeted genes", {
  expect_is(mytable, 'data.table')
  expect_equal(nrow(mytable), 2097)
  expect_equal(ncol(mytable), 9)
  expect_match(colnames(mytable)[1], "tx_name")
  expect_match(colnames(mytable)[9], "threeUTRs")
})

summary <- summarizeQueryRegions(queryRegions = queryRegions,
                                 txdbFeatures = features)
test_that("Summarizing query regions with txdbFeatures", {
  expect_is(summary, 'matrix')
  expect_match(colnames(summary), 'count')
  expect_equal(nrow(summary), length(names(features)))
})

profile <- calculateCoverageProfile(queryRegions = queryRegions,
                                   targetRegions = features[[1]],
                                   sampleN = 1000)
test_that("Getting coverage profile of query regions
          over a single gene feature from txdbFeatures", {
  expect_is(profile, 'data.frame')
  expect_equal(colnames(profile), c('coverage', 'bins'))
})

profileList <- calculateCoverageProfileList(queryRegions = queryRegions,
                                            targetRegionsList = features,
                                            sampleN = 1000)
test_that("Getting a list of coverage profiles
          over all gene features from txdbFeatures", {
  expect_is(profileList, 'list')
  expect_equal(length(names(profileList)), length(features))
})

profile <- calculateCoverageProfileFromTxdb(queryRegions = queryRegions,
                                                txdb = txdb,
                                                type = 'exons',
                                                sampleN = 1000)

test_that("Getting coverage profile of query regions
          over a single gene feature from a TxDb object", {
            expect_is(profile, 'data.frame')
            expect_equal(colnames(profile), c('coverage', 'bins'))
          })

profileList <- calculateCoverageProfileListFromTxdb(queryRegions = queryRegions,
                                                    txdb = txdb,
                                                    sampleN = 1000)
test_that("Getting a list of coverage profiles
          over all gene features from a txdb object", {
            expect_is(profileList, 'list')
            expect_equal(length(names(profileList)), length(features))
          })



