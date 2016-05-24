## ----global_options, include=FALSE---------------------------------------
knitr::opts_chunk$set(warning = FALSE, message = FALSE, eval = FALSE)

## ----load_libraries, results='hide'--------------------------------------
#  library(RCAS)

## ----sample_data---------------------------------------------------------
#  library(RCAS)
#  data(queryRegions) #sample queryRegions in BED format()
#  data(gff)          #sample GFF file

## ----RCAS_import_data----------------------------------------------------
#  #library(RCAS)
#  #queryRegions <- importBed(filePath = BED_file, sampleN = 10000)
#  #gff <- importGtf(filePath = GFF_file)

## ----queryGFF------------------------------------------------------------
#  overlaps <- queryGff(queryRegions = queryRegions, gff = gff)
#  #data.table is used to do quick summary operations
#  overlaps.dt <- data.table(as.data.frame(overlaps))

## ----query_gene_types----------------------------------------------------
#  biotype_col <- grep('biotype', colnames(overlaps.dt), value = TRUE)
#  df <- overlaps.dt[,length(unique(overlappingQuery)), by = biotype_col]
#  colnames(df) <- c("feature", "count")
#  df$percent <- round(df$count / length(queryRegions) * 100, 1)
#  df <- df[order(count, decreasing = TRUE)]
#  p <- plot_ly(data = df,
#               type = "bar",
#               x = feature,
#               y = percent,
#               text = paste("count:", count), color=feature)
#  layout(p = p,
#         margin = list(l=100, r=100, b=150),
#         xaxis = list(showticklabels = TRUE,  tickangle = 90),
#         yaxis = list(title = paste("percentage of query regions,",
#                                    "n =", length(queryRegions))))

## ----chromosomes_gene_features-------------------------------------------
#  df <- overlaps.dt[,length(unique(overlappingQuery)), by = c('seqnames', 'type')]
#  colnames(df) <- c('seqnames', 'type', 'count')
#  df <- df[order(seqnames)]
#  p <- plot_ly(data = df,
#               type = "bar",
#               x = seqnames,
#               y = count,
#               text = paste("count:", count),
#               color = type)
#  layout(p = p,
#         margin = list(b=150),
#         xaxis = list(showticklabels = TRUE,  tickangle = 90),
#         yaxis = list(title = paste("Number of query regions,",
#                                    "n =", length(queryRegions)
#                                    )
#                      )
#         )

## ----getTxdbFeatures-----------------------------------------------------
#  txdbFeatures <- getTxdbFeaturesFromGff(gff)

## ----summarizeQueryRegions-----------------------------------------------
#  summary <- summarizeQueryRegions(queryRegions = queryRegions,
#                                   txdbFeatures = txdbFeatures)
#  
#  df <- data.frame(summary)
#  df$percent <- round((df$count / length(queryRegions)), 3) * 100
#  p <- plot_ly( data = df,
#                x = rownames(df),
#                y = percent,
#                type = 'bar',
#                text = paste("count:", count),
#                color = rownames(df)
#                )
#  layout(p = p,
#         xaxis = list(title = 'features'),
#         yaxis = list(title = paste("percentage of query regions,",
#                                    "n =", length(queryRegions)
#                                    )
#                      ),
#         margin = list(b = 150, r = 50)
#         )

## ----getTargetedGenesTable-----------------------------------------------
#  dt = getTargetedGenesTable(queryRegions = queryRegions,
#                             txdbFeatures = txdbFeatures)
#  datatable(data = dt[order(transcripts, decreasing = TRUE)][1:10],
#            extensions = 'FixedColumns',
#            options = list(
#              dom = 't',
#              scrollX = TRUE,
#              scrollCollapse = TRUE
#              )
#            )

## ----coverageprofile-----------------------------------------------------
#  cov <- calculateCoverageProfile(queryRegions = queryRegions,
#                                 targetRegions = txdbFeatures$threeUTRs,
#                                 sampleN = 10000)
#  
#  p <- plot_ly(cov, x = bins, y = coverage)
#  p %>%   add_trace(y = fitted(loess(coverage ~ as.numeric(bins)))) %>%
#  layout(title = paste("Coverage along 3'UTRs (5' -> 3' direction)", sep = " "),
#         showlegend = FALSE,
#         margin = list(l = 50, r = 50, b = 50, t = 50))
#  

## ----coverageprofilelist-------------------------------------------------
#  covList <- calculateCoverageProfileList(queryRegions = queryRegions,
#                                         targetRegionsList = txdbFeatures,
#                                         sampleN = 10000)
#  
#  df <- do.call('cbind', covList)
#  df <- df[,!grepl(colnames(df), pattern = '*.bins')]
#  df$bins <- c(1:100)
#  colnames(df) <- gsub(pattern = ".coverage", replacement = "", x = colnames(df))
#  mdf <- reshape2::melt(df, id.vars = c('bins'))
#  colnames(mdf) <- c('bins', 'feature', 'coverage')
#  p = plot_ly(data = mdf, x = bins, y = coverage, color = feature)
#  layout(p)

## ----motif_analysis------------------------------------------------------
#  motifResults <- runMotifRG(queryRegions = queryRegions,
#                             genomeVersion = 'hg19',
#                             motifN = 2, nCores = 2)
#  
#  par(mfrow = c(1,2), mar = c(2,2,2,2))
#  for (i in 1:length(motifResults$motifs)) {
#    motifPattern <- motifResults$motifs[[i]]@pattern
#    motifRG::plotMotif(match = motifResults$motifs[[i]]@match$pattern,
#                       main = paste0('Motif-',i,': ',motifPattern),
#                       entropy = TRUE)
#  }

## ----motif_analysis_table------------------------------------------------
#  summary <- getMotifSummaryTable(motifResults)
#  datatable(summary, extensions = 'FixedColumns',
#            options = list(
#              dom = 't',
#              scrollX = TRUE,
#              scrollCollapse = TRUE
#              )
#            )

## ----GO analysis---------------------------------------------------------
#  
#  #get all genes from the GTF data
#  backgroundGenes <- unique(gff$gene_id)
#  #get genes that overlap query regions
#  targetedGenes <- unique(overlaps$gene_id)
#  
#  #run TopGO
#  goResults <- runTopGO(ontology = 'BP',
#                        species = 'human',
#                        backgroundGenes = backgroundGenes,
#                        targetedGenes = targetedGenes)
#  
#  datatable(data = goResults[1:10,],
#            extensions = 'FixedColumns',
#            options = list(dom = 't',
#                           scrollX = TRUE,
#                           scrollCollapse = TRUE
#                           )
#            )
#  

## ----msigdb_analysis-----------------------------------------------------
#  #msigDB <- parseMsigdb(msigdbFile)
#  data(msigDB)
#  msigdbResults <- runMSIGDB(msigDB = msigDB,
#      backgroundGenes = backgroundGenes, targetedGenes = targetedGenes)
#  
#  datatable(msigdbResults[1:10,],
#      extensions = 'FixedColumns',
#    options = list(
#      dom = 't',
#      scrollX = TRUE,
#      scrollCollapse = TRUE
#    ))
#  

