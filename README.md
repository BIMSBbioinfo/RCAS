# RCAS project

[![Build Status](https://travis-ci.org/BIMSBbioinfo/RCAS.svg?branch=master)](https://travis-ci.org/BIMSBbioinfo/RCAS)
![codecov.io](https://codecov.io/github/BIMSBbioinfo/RCAS/coverage.svg?branch=master)

## Introduction

RCAS is an R/Bioconductor package designed as a generic reporting tool for the
functional analysis of transcriptome-wide regions of interest detected by
high-throughput experiments. Such transcriptomic regions could be, for instance,
signal peaks detected by CLIP-Seq analysis for protein-RNA interaction sites,
RNA modification sites (alias the epitranscriptome), CAGE-tag locations, or any
other collection of query regions at the level of the transcriptome. RCAS
produces in-depth annotation summaries and coverage profiles based on the
distribution of the query regions with respect to transcript features (exons,
introns, 5’/3’ UTR regions, exon-intron boundaries, promoter regions). Moreover,
RCAS can carry out functional enrichment analyses and discriminative motif
discovery. RCAS supports all genome versions that are available in
`BSgenome::available.genomes`

## installation:

### Installing from [Bioconductor](http://bioconductor.org/packages/3.4/bioc/html/RCAS.html)
`if (!requireNamespace("BiocManager", quietly=TRUE))`
    `install.packages("BiocManager")`

`BiocManager::install('RCAS')`

### Installing the development version from Github
```
library('devtools')
devtools::install_github('BIMSBbioinfo/RCAS')
```

### Installing via Bioconda channel

`conda install bioconductor-rcas -c bioconda`

### Installing via Guix

`guix package -i r r-rcas`

## usage:

### Package vignettes and reference manual

For detailed instructions on how to use RCAS, please see: 
  - [package vignette for single sample analysis](http://bioconductor.org/packages/release/bioc/vignettes/RCAS/inst/doc/RCAS.vignette.html) 
  - [package vignette for multi-sample analysis](http://bioconductor.org/packages/release/bioc/vignettes/RCAS/inst/doc/RCAS.metaAnalysis.vignette.html)
  - [reference manual](http://bioconductor.org/packages/3.4/bioc/manuals/RCAS/man/RCAS.pdf) for more information about the detailed functions available in RCAS.

## Use cases from published RNA-based omics datasets

### Multi-sample analysis use case

- See an [example report](https://bimsbstatic.mdc-berlin.de/akalin/buyar/RCAS/1.4.2/RCAS.html) comparing the peak regions discovered via CLIP-sequencing experiments of the [RNA-binding protein FUS](http://www.uniprot.org/uniprot/P35637) by [Nakaya et al, 2013](https://www.ncbi.nlm.nih.gov/pubmed/23389473), [Synaptic Functional Regulator FMR1](http://www.uniprot.org/uniprot/Q06787) by [Ascano et al. 2012](https://www.ncbi.nlm.nih.gov/pubmed/23235829), and [Eukaryotic initiation factor 4A-III](http://www.uniprot.org/uniprot/P38919) by [Sauliere et al, 2012](https://www.ncbi.nlm.nih.gov/pubmed/23085716).


### Single Sample Analysis Use Cases

- [RCAS report](https://bimsbstatic.mdc-berlin.de/akalin/buyar/RCAS/1.1.1/PARCLIP_PUM2_Hafner2010b_hg19.bed.RCAS.report.html) for [PUM2](http://www.uniprot.org/uniprot/Q8TB72) RNA-binding sites detected by **PAR-CLIP** technique ([Hafner et al, 2010](https://www.ncbi.nlm.nih.gov/pubmed/20371350))
  - input:  [PARCLIP_PUM2_Hafner2010b_hg19](http://dorina.mdc-berlin.de/api/v1.0/download/regulator/hg19/PARCLIP_PUM2_Hafner2010b_hg19)


- [RCAS report](https://bimsbstatic.mdc-berlin.de/akalin/buyar/RCAS/1.1.1/PARCLIP_QKI_Hafner2010c_hg19.bed.RCAS.report.html) for [QKI](http://www.uniprot.org/uniprot/Q96PU8) RNA-binding sites detected by **PAR-CLIP** technique ([Hafner et al, 2010](https://www.ncbi.nlm.nih.gov/pubmed/20371350))
  - input: [PARCLIP_QKI_Hafner2010c_hg19](http://dorina.mdc-berlin.de/api/v1.0/download/regulator/hg19/PARCLIP_QKI_Hafner2010c_hg19 )


- [RCAS report](https://bimsbstatic.mdc-berlin.de/akalin/buyar/RCAS/1.1.1/PARCLIP_IGF2BP123_Hafner2010d_hg19.bed.RCAS.report.html) for IGF2BP[1](http://www.uniprot.org/uniprot/Q9NZI8)[2](http://www.uniprot.org/uniprot/Q9Y6M1)[3](http://www.uniprot.org/uniprot/O00425) RNA-binding sites detected by **PAR-CLIP** technique ([Hafner et al, 2010](https://www.ncbi.nlm.nih.gov/pubmed/20371350))
  - input:  [PARCLIP_IGF2BP123_Hafner2010d_hg19](http://dorina.mdc-berlin.de/api/v1.0/download/regulator/hg19/PARCLIP_IGF2BP123_Hafner2010d_hg19)


- [RCAS report](https://bimsbstatic.mdc-berlin.de/akalin/buyar/RCAS/1.1.1/human_FANTOM4_tiRNAs-hg19.bed.RCAS.report.html) for tiny RNA (tiRNA) loci detected by **deepCAGE** analysis ([Taft et al. 2009](https://www.ncbi.nlm.nih.gov/pubmed/19377478))
    - input: [human_FANTOM4_tiRNAs.bed](http://fantom.gsc.riken.jp/4/download/Supplemental_Materials/Taft_et_al_2009/human_FANTOM4_tiRNAs.bed)


- [RCAS report](https://bimsbstatic.mdc-berlin.de/akalin/buyar/RCAS/1.1.1/m1Asites.dominissini.2016.resized.bed.RCAS.report.html) for m<sup>1</sup>A **methylation sites** ([Dominissini et al, 2016](https://www.ncbi.nlm.nih.gov/pubmed/26863196))
    - input: [GSE70485_human_peaks.txt.gz](http://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE70485&format=file&file=GSE70485%5Fhuman%5Fpeaks%2Etxt%2Egz)

## Citation
In order to cite RCAS, please use: 

Bora Uyar, Dilmurat Yusuf, Ricardo Wurmus, Nikolaus Rajewsky, Uwe Ohler, Altuna Akalin; RCAS: an RNA centric annotation
  system for transcriptome-wide regions of interest. Nucleic Acids Res 2017 gkx120. doi: 10.1093/nar/gkx120

See our publication [here](https://academic.oup.com/nar/article/3038237/RCAS-an-RNA-centric-annotation-system-for). 

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
