require(knitr) # required for knitting from rmd to md
require(markdown) # required for md to html 
knit('./rcas.Rmd', 'rcas.md') # creates md file
markdownToHTML('rcas.md', 'rcas.html') # creates html file
#browseURL(paste('file://', file.path(getwd(),'rcas.html'), sep='')) # open file in browser 
