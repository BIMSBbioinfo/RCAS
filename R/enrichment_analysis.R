#' findEnrichedFunctions
#' 
#' Find enriched functional terms among the genes that overlap
#' the regions of interest. 
#' 
#' This function is basically a call to gprofiler2::gost function. 
#' It is here to serve as a replacement for other deprecated functional
#' enrichment functions. 
#' 
#' @examples 
#' data(gff)
#' data(queryRegions)
#' 
#' overlaps <- queryGff(queryRegions, gff)
#' res <- findEnrichedFunctions(unique(overlaps$gene_id), 'hsapiens')
#' 
#' @param targetGenes Vector of Ensembl gene ids or gene names 
#' @param species First letter of genus + species name: e.g. hsapiens
#' @param ... Other arguments to be passed to gprofiler2::gost 
#' @importFrom gprofiler2 gost
#' @export
findEnrichedFunctions <- function(targetGenes, species, ...) {
 res <- gprofiler2::gost(query = targetGenes, organism = species)
 return(res$result)
}

#' runTopGO
#' 
#' This function is deprecated. Use findEnrichedFunctions instead. 
#' @export
runTopGO <- function () {
  .Deprecated('findEnrichedFunctions')
}

#' parseMsigdb
#'
#' This function is deprecated. For functional enrichment 
#' analysis, use findEnrichedFunctions. 
#' 
#' @export
parseMsigdb <- function(){
  .Deprecated('findEnrichedFunctions.')
}

#' Print MSIGDB Dataset to a file
#' 
#' This function is deprecated. For functional enrichment 
#' analysis, use findEnrichedFunctions. 
#' @export
printMsigdbDataset = function(){
  .Deprecated('findEnrichedFunctions.')
}

#' retrieveOrthologs
#' 
#' This function is deprecated. For functional enrichment 
#' analysis, use findEnrichedFunctions. 
#'
#' @export
retrieveOrthologs <- function(){
  .Deprecated('findEnrichedFunctions')
}

#' createOrthologousMsigdbDataset
#' This function is deprecated. For functional enrichment 
#' analysis, use findEnrichedFunctions. 
#' @export
createOrthologousGeneSetList <- function() {
  .Deprecated('findEnrichedFunctions')
}

#' runGSEA
#' 
#' This function is deprecated. For functional enrichment 
#' analysis, use findEnrichedFunctions. 
#' 
#' @export
runGSEA <- function () {
  .Deprecated('findEnrichedFunctions')
}




