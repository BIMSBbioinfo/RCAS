library(RCAS)
context("tests for setting up and modifying sqlite database")

## prepare data for tests
if(file.exists('./test.RCASDB.sqlite')) {
  file.remove('./test.RCASDB.sqlite')
}

if(file.exists('./myProjDataFile.tsv')){
  file.remove('./myProjDataFile.tsv')
}

FUS_path <- system.file("extdata", "FUS_Nakaya2013c_hg19.bed",
                        package='RCAS')

FMR1_path <- system.file("extdata","FMR1_Ascano2012a_hg19.bed", 
                         package='RCAS')

projData <- data.frame('sampleName' = c('FUS', 'FMR1'), 
                       'bedFilePath' = c(FUS_path,FMR1_path), 
                       stringsAsFactors = FALSE)

write.table(projData, 
            'myProjDataFile.tsv', 
            sep = '\t', quote =FALSE,
            row.names = FALSE)

gtfFilePath <- system.file("extdata", "hg19.sample.gtf", package='RCAS')

# Testing createDB function
test_that("Testing createDB function", {
  expect_equal(object = createDB(dbPath = './test.RCASDB.sqlite', 
                        projDataFile = './myProjDataFile.tsv',
                        gtfFilePath = gtfFilePath,
                        genomeVersion = 'hg19',
                        motifAnalysis = FALSE,
                        coverageProfiles = FALSE), 
               expected = './test.RCASDB.sqlite')
  
  expect_error(createDB(dbPath = './test.RCASDB.sqlite', 
                        projDataFile = './myProjDataFile.tsv',
                        gtfFilePath = gtfFilePath,
                        genomeVersion = 'hg19', 
                        update = FALSE),
               regexp = 'A database already exists')
})

## Testing summarizeDatabaseContent function 
summary <- summarizeDatabaseContent(dbPath = './test.RCASDB.sqlite')
test_that("Testing summarizeDatabaseContent function", {
  expect_is(summary, 'data.frame')
  expect_equal(nrow(summary), 2)
  expect_equal(ncol(summary), 5)
  expect_identical(as.vector(summary$sampleName), c("FMR1", "FUS"))
})

mydb <- deleteSampleDataFromDB(dbPath = './test.RCASDB.sqlite', sampleNames = 'FUS')
samples <- RSQLite::dbReadTable(mydb, 'processedSamples')

## Testing deleteSampleDataFromDB function
test_that("Testing deleteSampleDataFromDB function", {
  expect_is(mydb, "SQLiteConnection")
  expect_false('FUS' %in% samples$sampleName)
  expect_true('FMR1' %in% samples$sampleName)
})

RSQLite::dbDisconnect(mydb)

if(file.exists('./test.RCASDB.sqlite')) {
  file.remove('./test.RCASDB.sqlite')
}

if(file.exists('./myProjDataFile.tsv')){
  file.remove('./myProjDataFile.tsv')
}
