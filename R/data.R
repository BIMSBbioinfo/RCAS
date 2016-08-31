#' Sample GFF file imported as a GRanges object
#'
#' This dataset contains genomic annotation data from Ensembl version 75 for
#' Homo sapiens downloaded from Ensembl. The GFF file is imported via the
#' \code{importGtf} function and a subset of the data is selected by choosing
#' features found on 'chr1'.
#'
#' @return A GRanges object
#' @format GRanges object with 238010 ranges and 16 metadata columns
#' @source \url{ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz}
"gff"

#' Sample BED file imported as a GRanges object
#'
#' This dataset contains a randomly selected sample of human LIN28A protein
#' binding sites detected by HITS-CLIP analysis downloaded from DoRina database
#' (LIN28A HITS-CLIP hESCs (Wilbert 2012)). The BED file is imported via the
#' \code{importBed} function and a subset of the data is selected by randomly
#' choosing 10000 regions.
#'
#' @return A GRanges object
#' @format GRanges object with 10000 ranges and 2 metadata columns
#' @source \url{http://dorina.mdc-berlin.de/regulators}
"queryRegions"


#' Random test gene sets
#'
#' This dataset contains random sets of genes with Entrez ids that is designed
#' to represent the data that can be parsed from MSIGDB database. 
#' using the \code{parseMsigdb} function. 
#' 
#' Actual curated datasets must be downloaded from the MSIGDB database
#'
#' @return A list object
#' @format A list of vectors, where each list element corresponds to a (randomized)
#' gene set, where genes are represented by Entrez ids.
"geneSets"
