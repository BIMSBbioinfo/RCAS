## ----load_libraries, results='hide', message=FALSE, warning=FALSE, echo=FALSE----
library(RCAS)
library(data.table)
library(plotly)
library(DT)
library(motifRG)

## ----importGtf, warning=FALSE, message=FALSE-----------------------------
#gff = importGtf(filePath = gtfFile)
gff <- readRDS('../data/hg19.chr1.gtf.granges.rds')

## ----importBed, warning=FALSE, message=FALSE-----------------------------
#queryRegions = importBed(filePath = peaksFile, sampleN = 10000)
queryRegions <- readRDS('../data/sample.bed.rds')

## ----queryGFF, warning=FALSE, message=FALSE------------------------------
overlaps = queryGff(queryRegions = queryRegions, gff = gff)
overlaps.dt = data.table(as.data.frame(overlaps)) #data.table is used to do quick summary operations 

## ----query_gene_types, warning=FALSE, message=FALSE----------------------
biotype_col = grep('biotype', colnames(overlaps.dt), value=T)
df = overlaps.dt[,length(unique(overlappingQuery)), by=biotype_col]
colnames(df) = c("feature", "count")
df$percent = round(df$count / length(queryRegions) * 100, 1)
df = df[order(count, decreasing = TRUE)]
p = plot_ly(df, type="bar", x = feature, y = percent, text=paste("count:", count), color=feature)
layout(p, margin = list(l=100, r=100, b=150), xaxis = list(showticklabels = TRUE,  tickangle = 90), yaxis = list(title = paste("percentage of query regions,", "n =",length(queryRegions))))

## ----chromosomes_gene_features, warning=FALSE, message=FALSE-------------
df = overlaps.dt[,length(unique(overlappingQuery)), by=c('seqnames', 'type')]
colnames(df) = c('seqnames', 'type', 'count')
df = df[order(seqnames)]
p = plot_ly(df, type="bar", x = seqnames, y = count, text=paste("count:", count), color=type)
layout(p, margin = list(b=150), xaxis = list(showticklabels = TRUE,  tickangle = 90), yaxis = list(title = paste("Number of query regions,", "n =",length(queryRegions))))

## ----getTxdbFeatures, warning=FALSE, message=FALSE-----------------------
txdbFeatures = getTxdbFeaturesFromGff(gff)

## ----getTargetedGenesTable, warning=FALSE, message=FALSE-----------------
dt = getTargetedGenesTable(queryRegions = queryRegions, txdbFeatures = txdbFeatures)
datatable(dt[order(transcripts, decreasing = T)][1:500], extensions = 'FixedColumns',
  options = list(
    dom = 't',
    scrollX = TRUE,
    scrollCollapse = TRUE
  ))

## ----coverageprofile, warning=FALSE, message=FALSE-----------------------
cov = calculateCoverageProfile(queryRegions = queryRegions, targetRegions = txdbFeatures$threeUTRs, sampleN = 10000)

p = plot_ly(cov, x = bins, y = coverage)
p %>%  
add_trace(y = fitted(loess(coverage ~ as.numeric(bins)))) %>%
layout(title = paste("Coverage along 3'UTRs (5' -> 3' direction)", sep=" "), showlegend = FALSE, margin = list(l= 50, r=50, b=50, t=50))


## ----coverageprofilelist, warning=FALSE, message=FALSE-------------------
covList = calculateCoverageProfileList(queryRegions = queryRegions, targetRegionsList = txdbFeatures, sampleN = 10000)

df = do.call('cbind', covList)
df = df[,!grepl(colnames(df), pattern = '*.bins')]
df$bins = c(1:100)
colnames(df) = gsub(pattern = ".coverage", replacement = "", x = colnames(df))
mdf = reshape2::melt(df, id.vars = c('bins'))
colnames(mdf) = c('bins', 'feature', 'coverage')
p = plot_ly(data = mdf, x=bins, y=coverage, color = feature)
layout(p)

## ----motif_analysis, warning=FALSE, message=FALSE------------------------
motifResults <- runMotifRG(queryRegions = queryRegions, genomeVersion = 'hg19', motifN = 4)

par(mfrow=c(2,2), mar=c(2,2,2,2))
for (i in 1:length(motifResults$motifs)){
  motifPattern = motifResults$motifs[[i]]@pattern
  motifRG::plotMotif(motifResults$motifs[[i]]@match$pattern, main=paste0('Motif-',i,': ',motifPattern), entropy = T)
}


## ----motif_analysis_table, warning=FALSE, message=FALSE------------------
summary <- getMotifSummaryTable(motifResults)
datatable(summary, extensions = 'FixedColumns',
  options = list(
    dom = 't',
    scrollX = TRUE,
    scrollCollapse = TRUE
  ))

## ----GO analysis, warning=FALSE, message=FALSE---------------------------

#get all genes from the gff data
backgroundGenes = unique(gff$gene_id)
#get genes that overlap query regions
targetedGenes = unique(overlaps$gene_id)

#run TopGO 
goResults = runTopGO(ontology = 'BP', species = 'human', backgroundGenes = backgroundGenes, targetedGenes = targetedGenes)

datatable(goResults[goResults$bh < 0.1,], extensions = 'FixedColumns',
  options = list(
    dom = 't',
    scrollX = TRUE,
    scrollCollapse = TRUE
  ))


## ----msigdb_analysis, warning=FALSE, message=FALSE-----------------------
#msigDB <- parseMsigdb(msigdbFile) 
msigDB <- readRDS('../data/msigdb.sample.rds')
msigdbResults <- runMSIGDB(msigDB = msigDB, backgroundGenes = backgroundGenes, targetedGenes = targetedGenes)

datatable(msigdbResults[msigdbResults$BH < 0.1,], extensions = 'FixedColumns',
  options = list(
    dom = 't',
    scrollX = TRUE,
    scrollCollapse = TRUE
  ))


