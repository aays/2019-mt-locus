# remove_duplicates.R
# 
# removes duplicates in the transposed fasta file
#
# usage:
# Rscript [infile.txt] [outfile.txt]

library(dplyr, warn.conflicts = FALSE)
library(readr)
library(magrittr)

args <- commandArgs(trailingOnly = TRUE)
infile <- args[1]
outfile <- args[2]

d <- read_delim(infile, delim = ' ', col_types = cols())
d %<>%
    distinct() %>%
    arrange(position) %>%
    mutate(previous_position = lag(position, default = -1)) %>%
    filter(position != previous_position) %>%
    select(-previous_position)

write_delim(d, path = outfile, delim = ' ')
