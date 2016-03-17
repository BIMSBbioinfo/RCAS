# RCAS -- standalone RNA centric annotation system
# Copyright © 2016 Ricardo Wurmus <rekado@elephly.net>
#
# This file is part of RCAS.
#
# RCAS is free software; see LICENSE file for details.

bin_SCRIPTS =					\
  src/RCAS

dist_libexec_SCRIPTS =				\
  src/rcas.Rmd					\
  src/rcas.msigdb.R				\
  src/rcas.GO.R					\
  src/rcas.motif.R				\
  src/run.rcas.R				\
  src/run.rcas.parallel.R			\
  src/map_msigdb_orthologs.R

basedir = $(datadir)/RCAS/base
dist_base_DATA = 				\
  src/base/c2.cp.v5.0.entrez.ce10.gmt		\
  src/base/c2.cp.v5.0.entrez.dm3.gmt		\
  src/base/c2.cp.v5.0.entrez.gmt		\
  src/base/c2.cp.v5.0.entrez.mm9.gmt		\
  src/base/header.html				\
  src/base/custom.css				\
  src/base/Caenorhabditis_elegans-U2T.meme	\
  src/base/Drosophila_melanogaster-U2T.meme	\
  src/base/Homo_sapiens-U2T.meme		\
  src/base/Mus_musculus-U2T.meme

baseimgdir = $(basedir)/img
dist_baseimg_DATA = src/base/img/**

dist_doc_DATA =					\
  README.md

# This is a hack to let the Python wrapper and thus snakemake know
# where the base files and all other libexec scripts are located.  All
# internal scripts called by snakemake rules and all snakefile parts
# are located in libexecdir.
src/RCAS: src/run.rcas.R.in
	sed -e "s|#:MAKE|\nlibexecdir = \"$(libexecdir)\"\nbasedir = \"$(basedir)\"|" -i $<
	mv $< $@

# Pandoc expects images to be relative to the current working
# directory when "--standalone" is provided.  Instead of hacking
# around this behaviour (by changing the working directory at runtime)
# we just patch in the full path to the basedir.

# Since this Makefile.am is processed by the configure script, we
# should not use "@basedir@" here, as that would be replaced with the
# actual path, leading to a very pointless "sed" command.
src/base/header.html: src/base/header.html.in
	sed -e "s|@basedir|$(basedir)|" $< > $@

CLEANFILES = src/base/header.html src/RCAS
