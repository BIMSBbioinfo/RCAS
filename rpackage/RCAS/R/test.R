#' @export
importGtf <- function (filePath, readFromRds = TRUE, saveObjectAsRds = TRUE, overwriteObjectAsRds = FALSE, keepStandardChr = TRUE) {

  rdsFilePath <- paste0(filePath, ".granges.rds")

  if (readFromRds == TRUE && file.exists(rdsFilePath)){
    cat('Reading existing granges.rds object from',rdsFilePath,'\n')
    gff <- readRDS(rdsFilePath)
  } else {
    cat('importing gtf file from', filePath,'\n')
    gff <- rtracklayer::import.gff(filePath)
  }

  if (!is.null(gff)) {
    GenomeInfoDb::seqlevelsStyle(gff) = 'UCSC'
    if (keepStandardChr == TRUE) {
      cat('Keeping standard chromosomes only\n')
      gff = GenomeInfoDb::keepStandardChromosomes(gff)
    }

    if (saveObjectAsRds == TRUE) {
      if (!file.exists(rdsFilePath)){
        cat('Saving gff object in RDS file at',rdsFilePath,'\n')
        saveRDS(gff, file=rdsFilePath)
      } else {
        if (overwriteObjectAsRds == TRUE) {
          cat('Overwriting gff object in RDS file at',rdsFilePath,'\n')
          saveRDS(gff, file=rdsFilePath)
        } else {
          message('File',rdsFilePath,'already exists. Use overwriteObjectAsRds = TRUE to overwrite the file\n')
        }
      }
    return (gff)
    }
  } else{
    stop("Couldn't import the gtf from",filePath,"\n")
  }
}







