# RCAS project

[![Build Status](https://travis-ci.org/BIMSBbioinfo/RCAS.svg?branch=master)](https://travis-ci.org/BIMSBbioinfo/RCAS)
[![codecov.io](https://codecov.io/github/BIMSBbioinfo/RCAS/coverage.svg?branch=master)]
Make a standalone RNA Centric Annotation System that
provides intuitive reports and publication ready graphics.

RCAS takes input peak intervals in BED foramt from clip-seq data
and automatically generates distributions of annotation features,
detected motifs, GO-term enrichment, pathway enrichment
and genomic coverage.

This repository provides the R package implementing the core features
of RCAS.  It can be used as a library.  For tools using the RCAS
library, such as a command line tool and a web interface, please check
the repository at <https://github.com/BIMSBbioinfo/RCAS-tools>.
