#' plotFeatureBoundaryCoverage
#' 
#' This function is used to create interactive plots displaying 5' and 3'
#' end coverage profiles of given transcript features. 
#' 
#' @param cvgF data.frame object containing 'fiveprime' coverage data returned 
#'   by getFeatureBoundaryCoverage function
#' @param cvgT data.fram object containing 'threeprime' coverage data returned 
#'   by getFeatureBoundaryCoverage function
#' @param featureName character object. This is used to label the axes (e.g.
#'   transcripts, exons)
#' @examples 
#' data(queryRegions)
#' data(gff)
#' txdb <- GenomicFeatures::makeTxDbFromGRanges(gff)
#' transcriptCoords <- GenomicFeatures::transcripts(txdb)
#' cvgF <- getFeatureBoundaryCoverage (queryRegions = queryRegions,
#'                                     featureCoords = transcriptCoords,
#'                                     flankSize = 100,
#'                                     boundaryType = 'fiveprime',
#'                                     sampleN = 1000)
#'                                     
#' cvgT <- getFeatureBoundaryCoverage (queryRegions = queryRegions,
#'                                     featureCoords = transcriptCoords,
#'                                     flankSize = 100,
#'                                     boundaryType = 'threeprime',
#'                                     sampleN = 1000)
#' p <- plotFeatureBoundaryCoverage(cvgF = cvgF, 
#'                                  cvgT = cvgT, 
#'                           featureName = 'transcript')   
#'                                                           
#' @return a plotly htmlwidget is returned
#' @export 
plotFeatureBoundaryCoverage <- function (cvgF, cvgT, featureName) {

  p1 <- plotly::plot_ly(data = cvgF, x = ~bases, y = ~meanCoverage, 
                        type = 'scatter', mode = 'lines', 
                        name = paste(featureName, "5' end coverage")) %>% 
    add_ribbons(x = ~bases, 
                ymin = cvgF$meanCoverage - cvgF$standardError*1.96, 
                ymax = cvgF$meanCoverage + cvgF$standardError*1.96,
                line = list(color = 'rgba(182, 7, 127, 0.05)'),
                fillcolor = 'rgba(182, 7, 127, 0.22)',
                name = paste(featureName, "5' standard error (95% conf. int.)")) 
  
  p2 <- plotly::plot_ly(data = cvgT, x = ~bases, y = ~meanCoverage, 
                        type = 'scatter', mode = 'lines', 
                        name = paste(featureName, "3' end coverage")) %>% 
    add_ribbons(x = ~bases, 
                ymin = cvgT$meanCoverage - cvgT$standardError*1.96, 
                ymax = cvgT$meanCoverage + cvgT$standardError*1.96, 
                line = list(color = 'rgba(7, 164, 181, 0.05)'),
                fillcolor = 'rgba(7, 164, 181, 0.2)',
                name = paste(featureName, "3' standard error (95% conf. int.)")
    )
  
  p <- plotly::subplot(p1, p2, shareY = TRUE) %>% 
    layout (xaxis = list(title = "Distance (bp) to 5' boundary"), 
            xaxis2 = list(title = "Distance (bp) to 3' boundary"),
            yaxis = list(title = "Mean Coverage Score"),
            legend = list(x = 0, y = 100, orientation = 'h'),
            font = list(size = 14), 
            margin = list(l = 100, b = 100))
  return(p)
}