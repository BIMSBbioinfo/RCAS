# RCAS project

[![Build Status](https://travis-ci.org/BIMSBbioinfo/RCAS.svg?branch=master)](https://travis-ci.org/BIMSBbioinfo/RCAS)
![codecov.io](https://codecov.io/github/BIMSBbioinfo/RCAS/coverage.svg?branch=master)

RCAS is an R/Bioconductor package designed as a generic reporting tool
for the functional analysis of transcriptome-wide regions of interest detected
by high-throughput experiments. Such transcriptomic regions could be,
for instance, signal peaks detected by CLIP-Seq analysis for protein-RNA
interaction sites, RNA modification sites (alias the epitranscriptome),
CAGE-tag locations, or any other collection of query regions at the level of
the transcriptome. RCAS produces in-depth annotation summaries and
coverage profiles based on the distribution of the query regions with respect
to transcript features (exons, introns, 5’/3’ UTR regions, exon-intron
boundaries, promoter regions). Moreover, RCAS can carry out functional
enrichment analyses of annotated gene sets, GO terms, and de novo motif
discovery. RCAS is available in the Bioconductor repository, packaged in multiple
environments including Conda, Galaxy, and Guix, and as a webservice
at http://rcas.mdc-berlin.de/.

Currently supported genome builds are hg19 (human), mm9 (mouse), dm3 (fly), and ce10 (worm). 
Modules for annotation summaries and motif analysis are supported for each of these genome  builds. 
GO term and gene-set enrichment analyses are supported for hg19, mm9, and dm3 (ce10 is currently not supported 
for GO/GSEA modules). 

## installation:

### Installing from [Bioconductor](http://bioconductor.org/packages/3.4/bioc/html/RCAS.html)
`source('http://bioconductor.org/biocLite.R')`

`biocLite('RCAS')`

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

- GTF file (Required)

Firstly, download the input datasets from ENSEMBL (GTF file) and doRiNA (BED file), and optionally gene-set annotations from the [MSIGDB](http://software.broadinstitute.org/gsea/msigdb/)

```
wget ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz
gunzip Homo_sapiens.GRCh37.75.gtf.gz
```

- BED file (Required)

Sample BED file from a HITSCLIP study of LIN28A protein [Wilbert ML, et al, 2012](https://www.ncbi.nlm.nih.gov/pubmed/22959275) downloaded from DoRiNA database
```
wget https://bimsbstatic.mdc-berlin.de/akalin/buyar/RCAS/HITSCLIP_LIN28AWilbert2012a_hg19.bed
```

- Gene set annotations (Optional)

Download the human gene set annotations from the Molecular Signatures Database.
The recommended set of gene set annotations is the [curated gene sets](http://software.broadinstitute.org/gsea/msigdb/collections.jsp#C2) with 'Entrez' gene identifiers.

If the gene set enrichment analysis using the MSIGDB annotations is needed for other species, RCAS can take the MSIGDB annotations for human and map them to mouse and fly genomes via orthologous genes. In any case, it must be the human gene-set annotations provided to the main reporting function of RCAS (**runReport**). The mapping to other species happens in the background.

Generate an annotation report using the GTF file and built-in test data (BED format) as input.
```
library('RCAS')
RCAS::runReport(gffFilePath = "./Homo_sapiens.GRCh37.75.gtf",  queryFilePath = "./HITSCLIP_LIN28AWilbert2012a_hg19.bed", msigdbFilePath = <path to MSIGDB gene-sets with ENTREZ gene ids>, genomeVersion = 'hg19')
```

## Package vignette and reference manual

See the [package vignette](http://bioconductor.org/packages/3.4/bioc/vignettes/RCAS/inst/doc/RCAS.vignette.html) and [the reference manual](http://bioconductor.org/packages/3.4/bioc/manuals/RCAS/man/RCAS.pdf) for more information about the detailed functions available in RCAS.

## Use cases from published RNA-based omics datasets

- The reference files needed for reproducing the reports below:
    - [Homo_sapiens.GRCh37.75.gtf.gz]( http://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz)
    - [c2.all.v5.1.entrez.gmt](http://software.broadinstitute.org/gsea/msigdb/download_file.jsp?filePath=/resources/msigdb/5.1/c2.all.v5.1.entrez.gmt)


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
