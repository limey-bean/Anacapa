#!/usr/bin/env Rscript
# Interpret the input variables -----


args = commandArgs(trailingOnly=TRUE)
# check that there are 3 args
if (length(args) != 3) {
  stop("please make sure there are 3 arguments: the path to the text file to be summarized, the barcode name, the path to the output file")
}


# need base library and dplyr, if ran dada2 should have it installed
# but just incase you do not have dplyr... Download package from CRAN
.cran_packages  <-  c("dplyr","readr", "stringr")
.inst <- .cran_packages %in% installed.packages()
if (any(!.inst)) {
  install.packages(.cran_packages[!.inst], repos = "http://cran.rstudio.com/")
}

inputfile = args[1]
metaB = args[2]
outputfile = args[3]



library("dplyr")
library("readr")
library("stringr")

# open file to summarize
df <- read_delim(inputfile, delim = "\t" )

# ignore the first column and for each row with the same taxonomic path, sum the counts per column
by_taxon <- df %>% mutate(sum.taxonomy = stringr::str_replace_all(as.character(`sum.taxonomy`), ";;", ";NA;")) %>%
  group_by(sum.taxonomy) %>% summarize_if(is.numeric, sum)

# export summarized file

write.table(by_taxon, file = outputfile, row.names=FALSE, sep="\t", quote=FALSE)
