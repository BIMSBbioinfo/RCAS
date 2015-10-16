
#RCAS project
Make a standalone RNA Centric Annotation System that provides intuitive reports and publication ready graphics.

##Tutorial

###RCAS dependencies:

-**snakemake**: https://bitbucket.org/johanneskoester/snakemake/wiki/Documentation

-**bedtools**:
http://bedtools.readthedocs.org/en/latest/content/installation.html

-**parse_anot.py**: src/parse_anot.py

-**gff track** e.g.  wget -qO- ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_23/gencode.v23.annotation.gff3.gz     | gunzip > gencode.v23.annotation.gff3

###RCAS workflow is built with snakemake:

-rules: src/RCAS.snakefile

-visualization: RCAS_dag.pdf

###Test case
####Aim: **generate annotation** for coordinates in the target bed files

1. clone the source: **git clone https://github.com/BIMSBbioinfo/RCAS**

2. make dependency tools system wise executable, e.g. **ln -s path_to/RCAS/src/parse_anot.py ~/bin/.**

3. cd path_to/RCAS/test

4. tartget files in path_to/RCAS/test: **PARCLIP_AGO1234_Hafner2010a_hg19_xaa.bed  PARCLIP_AGO1234_Hafner2010a_hg19_xab.bed  PARCLIP_AGO1234_Hafner2010a_hg19_xac.bed** 

5. command: **snakemake -s ../src/RCAS.snakefile -p --config ref=path_to/gencode.v23.annotation.gff3 infile=PAR*bed**

6. output: **PARCLIP_AGO1234_Hafner2010a_hg19_xaa.anot.tsv  PARCLIP_AGO1234_Hafner2010a_hg19_xab.anot.tsv  PARCLIP_AGO1234_Hafner2010a_hg19_xac.anot.tsv**





    
