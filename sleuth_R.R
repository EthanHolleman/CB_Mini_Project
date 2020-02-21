#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)
library(sleuth)
library(dplyr)

stab <- read.csv(args[0], header=TRUE, stringsAsFactors = FALSE)
so <- sleuth_prep(stab)
so <- sleuth_fit(so, ~condition, 'full')

# extract sleuth results

sleuth_t <- sleuth_results(so, 'reduced:full', 'Irt', show_all=FALSE)
sig_results <- dplyr::filter(sleuth_t, qval <= 0.05) %>%dplyr::arrange(pval)

# write table of significant results to a new csv file
# using the second arguement as the write path

write.table(sig_results, file=args[1], quote = FALSE, row.names = FALSE)