#RCAS project
Make a standalone RNA Centric Annotation System that
provides intuitive reports and publication ready graphics.

RCAS takes input peak intervals in BED foramt from clip-seq data
and automatically generates distributions of annotation features,
detected motifs, GO-term enrichment, pathway enrichment
and genomic coverage.

The workflow is modularized so user can switch on/off
optional steps including motif detection, GO-term enrichment,
pathway enrichment and genomic coverage.
By default these optional steps are switched off
and RCAS only generates distributions of annotation features
including biotypes and genomic features.

###RCAS dependencies:

-**snakemake**: https://bitbucket.org/johanneskoester/snakemake/wiki/Documentation

-**bedtools**:
http://bedtools.readthedocs.org/en/latest/content/installation.html

--**fastaFromBed** (from bedtools)

-**MEME-chip**:
http://meme-suite.org/meme-software/4.10.2/meme_4.10.2.tar.gz

-**pandoc (>= 1.12.3)**

-**rmarkdown**

-**rtracklayer**

-**data.table**

-**biomaRt**

-**org.Hs.eg.db**

-**org.Ce.eg.db**

-**org.Mm.eg.db**

-**org.Dm.eg.db**

-**topGO**

-**DT**

-**plotly**

-**dplyr**

-**genomation**

-**GenomicFeatures**

-**gff track** e.g.  wget -qO- ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gff3.gz | gunzip > gencode.v19.annotation.gff3

-**genome reference** e.g. http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz

-**c2.cp.v5.0.entrez.gmt**:
http://software.broadinstitute.org/gsea/msigdb/download_file.jsp?filePath=/resources/msigdb/5.0/c2.cp.v5.0.entrez.gmt

###RCAS workflow is built with snakemake:

-rules of snakemake in src:

    RCAS.snakefile

-RCAS scripts in src:

    RCAS.py

    parse_anot.py

    top_motifs.py

    rcas.GO.R

    rcas.msigdb.R

    rcas.Rmd

    make.R

-Other dependencies in base:

    Homo_sapiens-U2T.meme

    Mus_musculus-U2T.meme

    Caenorhabditis_elegans-U2T.meme

    Drosophila_melanogaster-U2T.meme

    c2.cp.v5.0.entrez.gmt

    c2.cp.v5.0.entrez.dm3.gmt

    c2.cp.v5.0.entrez.mm9.gmt

    c2.cp.v5.0.entrez.ce10.gmt

    custom.css

    header.html

    img folder

###Test case
####Aim: generate analysis report for intervals in the target BED files

*Besides the default annotation step, pathway enrichment is also enabled vi '--run_PATHrich True'*

1. clone the source: **git clone -b modular https://github.com/BIMSBbioinfo/RCAS**

2. cd /path/to/RCAS/test

3. tartget files in /path/to/RCAS/test: **xaa.bed  xab.bed**

4. check help message: **python2 ../src/RCAS.py -h**

4. sample command: **python2 ../src/RCAS.py --genome /path/to/hg19.fa --gff3 /path/to/gencode.v23.annotation.gff3 --RCAS_path ../ ../test/xaa.bed ../test/xab.bed --run_PATHrich True**

5. output: **xaa.rcas.html xab.rcas.html**
