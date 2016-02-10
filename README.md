# RCAS project

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

## Installing and building from source

RCAS uses autotools to provide a standard configuration script and
installation rules following the GNU coding standards.  To build from
a simple repository checkout you need to first bootstrap the
configuration script.

~~~
git clone https://github.com/BIMSBbioinfo/RCAS.git
cd RCAS
autoreconf -vif
./configure --prefix=/where/you/want/rcas
make install
~~~

When building from a release tarball this bootstrap step is not
required:

~~~
tar xf rcas-0.0.1.tar.gz
cd rcas
./configure --prefix=/where/you/want/rcas
make install
~~~

If your system provides both Python 2.x and Python 3.x and defaults to
Python 3.x, you may need to pass `PYTHON=/path/to/python2` to the
`configure` script to ensure that the Python 2.x interpreter is used,
for example:

~~~
./configure PYTHON=/bin/python2 --prefix=/opt/rcas
sudo make install
~~~

After the first official release we will also provide a Guix package
definition for RCAS so that you can install the pipeline with a single
command such as:

    guix package -i rcas


## Hacking on RCAS

For development purposes we also provide a package expression for GNU
Guix.  To enter a development shell in which all dependent tools are
available run this:

    guix environment -l package.scm

In this shell you can configure and install RCAS as usual, but all
tools and R packages will be handled by Guix.  There is no need to
manually install the tools.


## RCAS dependencies

We recommend to install dependencies with a package manager using one
of the package recipes we provide (see "Hacking on RCAS" above).  If
you really *must* install the dependencies manually please see the
following list for pointers:

### Tools

- [snakemake](https://bitbucket.org/johanneskoester/snakemake/wiki/Documentation)
- [R](https://www.r-project.org/)
- fastaFromBed from [bedtools](http://bedtools.readthedocs.org/en/latest/content/installation.html)
- [MEME-chip](http://meme-suite.org/meme-software/4.10.2/meme_4.10.2.tar.gz)
- pandoc (>= 1.12.3)

### R packages

- rmarkdown
- rtracklayer
- data.table
- biomaRt
- org.Hs.eg.db
- org.Ce.eg.db
- org.Mm.eg.db
- org.Dm.eg.db
- topGO
- DT
- plotly
- dplyr
- genomation
- GenomicFeatures

### Data files

- gff track e.g.  `wget -qO- ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gff3.gz | gunzip > gencode.v19.annotation.gff3`
- genome reference e.g. this
  [hg19 reference](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz)
- [c2.cp.v5.0.entrez.gmt](http://software.broadinstitute.org/gsea/msigdb/download_file.jsp?filePath=/resources/msigdb/5.0/c2.cp.v5.0.entrez.gmt)

## RCAS workflow is built with snakemake:

Note that files containing placeholders that are substituted at
configuration time end on `.in`.  Successful configuration generates
the target files without the `.in` extension.  In the lists below we
omit `.in` for clarity.

- The pipeline wrapper and entry point:
  - `src/RCAS`

- All snakemake rules end on `.snakefile`.
  - `src/RCAS.snakefile`
  - `src/RCAS.anot.snakefile`
  - `src/RCAS.motif.snakefile`
  - `src/RCAS.GOrich.snakefile`
  - `src/RCAS.PATHrich.snakefile`

- Custom scripts that are called at runtime:
  - `src/parse_anot.py`
  - `src/top_motifs.py`
  - `src/rcas.GO.R`
  - `src/rcas.msigdb.R`
  - `src/rcas.Rmd`
  - `src/generate_report.sh`

- Plain data files are located in `src/base`:
  - `Homo_sapiens-U2T.meme`
  - `Mus_musculus-U2T.meme`
  - `Caenorhabditis_elegans-U2T.meme`
  - `Drosophila_melanogaster-U2T.meme`
  - `c2.cp.v5.0.entrez.gmt`
  - `c2.cp.v5.0.entrez.dm3.gmt`
  - `c2.cp.v5.0.entrez.mm9.gmt`
  - `c2.cp.v5.0.entrez.ce10.gmt`
  - `custom.css`
  - `header.html`
  - `img` folder

## Test case

Aim: generate analysis report for intervals in the target BED files.

*Besides the default annotation step, pathway enrichment is also enabled vi '--run_PATHrich'*

1. Install as per the instructions above, e.g.
    ~~~
    autoreconf -vif
    ./configure --prefix=/opt/rcas
    make clean
    sudo make install
    ~~~

2. `cd /opt/rcas`

3. target files in `./test`: `xaa.bed`, `xab.bed`

4. check help message: `./bin/RCAS -h`

5. sample command:
    ~~~
    ./bin/RCAS --genome /path/to/hg19.fa                   \
               --gff3 /path/to/gencode.v19.annotation.gff3 \
               --run_PATHrich                              \
               ./test/xaa.bed ./test/xab.bed
    ~~~

6. output: `xaa.rcas.html`, `xab.rcas.html`
