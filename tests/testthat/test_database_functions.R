library(RCAS)
context("tests for setting up and modifying sqlite database")

dbPath <- system.file("extdata", "hg19.RCASDB.sqlite", package='RCAS')
dbPathTest <- file.path(getwd(), 'test.RCASDB.sqlite')
file.copy(from = dbPath, to = dbPathTest, overwrite = TRUE)

# Testing createDB function
test_that("Testing createDB function", {
  expect_error(createDB(dbPath = dbPath,
                        genomeVersion = 'hg19',
                        update = FALSE),
               regexp = 'A database already exists')
})

## Testing summarizeDatabaseContent function
summary <- summarizeDatabaseContent(dbPath = dbPathTest)
test_that("Testing summarizeDatabaseContent function", {
  expect_is(summary, 'data.frame')
  expect_equal(nrow(summary), 2)
  expect_equal(ncol(summary), 5)
  expect_identical(as.vector(summary$sampleName), c("FMR1", "FUS"))
})

mydb <- deleteSampleDataFromDB(dbPath = dbPathTest, sampleNames = 'FUS')
samples <- RSQLite::dbReadTable(mydb, 'processedSamples')

## Testing deleteSampleDataFromDB function
test_that("Testing deleteSampleDataFromDB function", {
  expect_is(mydb, "SQLiteConnection")
  expect_false('FUS' %in% samples$sampleName)
  expect_true('FMR1' %in% samples$sampleName)
})

RSQLite::dbDisconnect(mydb)

if(file.exists(dbPathTest)) {
  file.remove(dbPathTest)
}
