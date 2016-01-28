
#RCAS project
Make a standalone RNA Centric Annotation System that provides intuitive reports and publication ready graphics.

##Tutorial

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

-In-house scripts in src:

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
####Aim: **generate analysis report** for coordinates in the target bed files

1. clone the source: **git clone https://github.com/BIMSBbioinfo/RCAS**

2. cd path_to/RCAS/test

3. tartget files in path_to/RCAS/test: **PARCLIP_AGO1234_Hafner2010a_hg19_xaa.bed  PARCLIP_AGO1234_Hafner2010a_hg19_xab.bed  PARCLIP_AGO1234_Hafner2010a_hg19_xac.bed** 

4. sample command: **snakemake -s ../src/RCAS.snakefile -p --config gff3=path_to/gencode.v19.annotation.gff3 genome=path_to/hg19.fa RCAS_path=../ infile=PARCLIP_AGO1234_Hafner2010a_hg19_xaa.bed**

5. output: **PARCLIP_AGO1234_Hafner2010a_hg19_xaa.rcas.html**

