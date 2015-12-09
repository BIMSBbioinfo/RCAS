library('rmarkdown')

args <- commandArgs(TRUE)## Parse arguments (we expect the form --arg=value)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
argsDF <- as.data.frame(do.call("rbind", parseArgs(args)))
argsL <- as.list(as.character(argsDF$V2))
names(argsL) <- argsDF$V1

out_file = argsL$out
rmarkdown::render('rcas.Rmd', output_file = out_file)
