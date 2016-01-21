#!/bin/bash
echo "The script you are running has basename `basename $0`, dirname `dirname $0`"
echo "The present working directory is `pwd`"

SCRIPT_DIR=`dirname $0`
report_script="${SCRIPT_DIR}/rcas.Rmd"
css="${SCRIPT_DIR}/base/custom.css"
outdir=`pwd`
header="${SCRIPT_DIR}/base/header.html"

output_filename=$1

Rscript -e "library('rmarkdown'); rmarkdown::render('${report_script}', output_file = '${output_filename}', output_dir='${outdir}', html_document(toc=TRUE, theme='cerulean', number_sections=TRUE, css='${css}', includes=includes(before_body='${header}')))" ${outdir} ${2} ${3} ${4} ${5} ${6} ${7} ${8} ${9} ${10}
