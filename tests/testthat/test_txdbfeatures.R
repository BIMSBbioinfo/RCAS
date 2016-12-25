library(RCAS)
context("Functions to get txdbFeatures and calculating coverage profiles")

data(gff)
txdb <- GenomicFeatures::makeTxDbFromGRanges(gff)
test_that("Importing txdb", {
  expect_is(txdb, 'TxDb')
})

features <- getTxdbFeatures(txdb)
test_that("Parsing features from txdb", {
  expect_is(features, 'list')
  expect_equal(length(names(features)), 7)
  expect_match(names(features)[1], "transcripts")
  expect_match(names(features)[7], "threeUTRs")
})

features <- getTxdbFeaturesFromGRanges(gffData = gff)
test_that("Parsing features from GFF", {
  expect_is(features, 'list')
  expect_equal(length(names(features)), 7)
  expect_match(names(features)[1], "transcripts")
  expect_match(names(features)[7], "threeUTRs")
})

data(queryRegions)
mytable <- getTargetedGenesTable(queryRegions = queryRegions,
                                 txdbFeatures = features)
test_that("Getting table of targeted genes", {
  expect_is(mytable, 'data.table')
  expect_equal(nrow(mytable), 2094)
  expect_equal(ncol(mytable), 8)
  expect_match(colnames(mytable)[1], "tx_name")
  expect_match(colnames(mytable)[8], "threeUTRs")
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
  expect_is(profile, 'ScoreMatrix')
  expect_equal(ncol(profile), 100)
})

profileList <- calculateCoverageProfileList(queryRegions = queryRegions,
                                            targetRegionsList = features,
                                            sampleN = 1000)
test_that("Getting a list of coverage profiles
          over all gene features from txdbFeatures", {
  expect_is(profileList, 'data.frame')
  expect_equal(colnames(profileList), c('bins', 'meanCoverage', 'standardError', 'feature'))
  expect_equal(dim(profileList), c(700, 4))
})

profile <- calculateCoverageProfileFromTxdb(queryRegions = queryRegions,
                                                txdb = txdb,
                                                type = 'exons',
                                                sampleN = 1000)

test_that("Getting coverage profile of query regions
          over a single gene feature from a TxDb object", {
            expect_is(profile, 'ScoreMatrix')
            expect_equal(ncol(profile), 100)
          })

profileList <- calculateCoverageProfileListFromTxdb(queryRegions = queryRegions,
                                                    txdb = txdb,
                                                    sampleN = 1000)
test_that("Getting a list of coverage profiles
          over all gene features from a txdb object", {
            expect_is(profileList, 'list')
            expect_equal(length(names(profileList)), length(features))
          })

transcriptEndCoverage <- getFeatureBoundaryCoverage(
                                queryRegions = queryRegions,
                                featureCoords = features$transcripts[1:1000],
                                flankSize = 100, 
                                boundaryType = 'threeprime',
                                sampleN = 0
                                )

test_that("Calculating per base coverage profile at TES boundary", {
  expect_equal(colnames(transcriptEndCoverage), c("bases","meanCoverage","standardError"))
  expect_is(transcriptEndCoverage, 'data.frame')
  expect_equal(sum(transcriptEndCoverage$bases), -100)
})

transcriptEndCoverageBin <- getFeatureBoundaryCoverageBin(
  queryRegions = queryRegions,
  featureCoords = features$transcripts[1:1000],
  flankSize = 100,
  sampleN = 0
)

test_that("Calculating per bin coverage profile at TSS and TES boundaries", {
  expect_equal(colnames(transcriptEndCoverageBin), c("fivePrime","threePrime","bins"))
  expect_is(transcriptEndCoverageBin, 'data.frame')
  expect_equal(as.vector(colSums(transcriptEndCoverageBin)), c(15, 241, 0))
})


