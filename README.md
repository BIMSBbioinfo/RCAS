# RCAS project

[![Build Status](https://travis-ci.org/BIMSBbioinfo/RCAS.svg?branch=master)](https://travis-ci.org/BIMSBbioinfo/RCAS)
![codecov.io](https://codecov.io/github/BIMSBbioinfo/RCAS/coverage.svg?branch=master)

Make a standalone RNA Centric Annotation System that
provides intuitive reports and publication ready graphics.

RCAS takes as input transcript-level intervals of interest (for instance CLIP-Seq peaks, epitranscriptome modification sites, 
CAGE tags) in BED format and a GTF file for genomic feature annotations,
and automatically generates distributions of annotation features,
detected motifs, GO-term enrichment, gene-set enrichment and coverage profiles of signals across gene features such as exons, introns,
promoters, UTRs, exon-intron boundaries. The main function **runReport** can generate an HTML report with interactive figures and tables. 

This repository provides the R package implementing the core features
of RCAS.  It can be used as a library.  For tools using the RCAS
library, such as a command line tool and a web interface, please check
the repository at <https://github.com/BIMSBbioinfo/RCAS-tools>.

## installation:

### Installing from [Bioconductor](http://bioconductor.org/packages/3.4/bioc/html/RCAS.html) 
`source('http://bioconductor.org/biocLite.org')`

`biocLite('RCAS')`

### Installing the development version from Github
```
library('devtools')
devtools::install_github('BIMSBbioinfo/RCAS')
```

## usage: 

Firstly, download the input datasets from ENSEMBL (GTF file) and doRiNA (BED file), and optionally gene-set annotations from the [MSIGDB](http://software.broadinstitute.org/gsea/msigdb/)

- GTF file: 
```
wget ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz
gunzip Homo_sapiens.GRCh37.75.gtf.gz
```

- BED file
Sample BED file from a HITSCLIP study of LIN28A protein [Wilbert ML, et al, 2012](https://www.ncbi.nlm.nih.gov/pubmed/22959275) downloaded from DoRiNA database
```
wget https://bimsbstatic.mdc-berlin.de/akalin/buyar/RCAS/HITSCLIP_LIN28AWilbert2012a_hg19.bed
```

- Gene set annotations 
Optionally; download the human gene set annotations from the Molecular Signatures Database. 
The recommended set of gene set annotations is the [curated gene sets](http://software.broadinstitute.org/gsea/msigdb/collections.jsp#C2) with 'Entrez' gene identifiers. 

If the gene set enrichment analysis using the MSIGDB annotations is needed for other species, RCAS can take the MSIGDB annotations for human and map them to mouse and fly genomes via orthologous genes. In any case, it must be the human gene-set annotations provided to the main reporting function of RCAS (**runReport**). The mapping to other species happens in the background. 

Generate an annotation report using the GTF file and built-in test data (BED format) as input. 
```
library('RCAS')
RCAS::runReport(gffFilePath = "./Homo_sapiens.GRCh37.75.gtf",  queryFilePath = "./HITSCLIP_LIN28AWilbert2012a_hg19.bed", msigdbFilePath = <path to MSIGDB gene-sets with ENTREZ gene ids>, genomeVersion = 'hg19')
```

## Package vignette and reference manual

See the [package vignette](http://bioconductor.org/packages/3.4/bioc/vignettes/RCAS/inst/doc/RCAS.vignette.html) and [the reference manual](http://bioconductor.org/packages/3.4/bioc/manuals/RCAS/man/RCAS.pdf) for more information about the detailed functions available in RCAS. 

## Acknowledgements

RCAS is developed in the group of
[Altuna Akalin](http://bioinformatics.mdc-berlin.de/team.html#altuna-akalin-phd)
(head of the Scientific Bioinformatics Platform) by
[Bora Uyar](http://bioinformatics.mdc-berlin.de/team.html#bora-uyar-phd)
(Bioinformatics Scientist),
[Dilmurat Yusuf](http://bioinformatics.mdc-berlin.de/team.html#dilmurat-yusuf-phd)
(Bioinformatics Scientist) and
[Ricardo Wurmus](http://bioinformatics.mdc-berlin.de/team.html#ricardo-wurmus)
(System Administrator) at the Berlin Institute of Medical Systems Biology
([BIMSB](https://www.mdc-berlin.de/13800178/en/bimsb))
at the Max-Delbrueck-Center for Molecular Medicine
([MDC](https://www.mdc-berlin.de)) in Berlin.

RCAS is developed as a bioinformatics service as part of
the [RNA Bioinformatics Center](http://www.denbi.de/index.php/rbc),
which is one of the eight centers of
the German Network for Bioinformatics Infrastructure
([de.NBI](http://www.denbi.de/)).
