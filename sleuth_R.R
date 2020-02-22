#!/usr/bin/env Rscript
library("optparse")
library(sleuth)
library(dplyr)

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="out.txt", 
              help="output file name [default= %default]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

stab <- read.csv(opt$file, header=TRUE, stringsAsFactors = FALSE)
so <- sleuth_prep(stab)
so <- sleuth_fit(so, ~condition, 'full')
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')

# extract sleuth results

sleuth_t <- sleuth_results(so, 'reduced:full', 'lrt', show_all=FALSE)
sig_results <- dplyr::filter(sleuth_t, qval <= 0.05) %>%dplyr::arrange(pval)

# write table of significant results to a new csv file
# using the second arguement as the write path

write.table(sig_results, file=opt$out, quote = FALSE, row.names = FALSE)