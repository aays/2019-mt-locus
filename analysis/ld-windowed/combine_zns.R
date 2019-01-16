# combine_zns.R - combine zns files
#
# usage:
# Rscript script.R --dir [directory] --outfile [file]

library(readr)
library(dplyr, warn.conflicts = FALSE)
library(purrr)
library(magrittr, warn.conflicts = FALSE)
library(fs)
library(optparse)

args <- function() {
    option_list <- list(
        make_option(c("-d", "--directory"), type = "character", default = NULL, 
                  help = "Directory of ZnS files", metavar = "character"),
        make_option(c("-o", "--outfile"), type = "character", default = NULL, 
                  help= "File to write to", metavar = "character")
    ) 
    opt_parser <- OptionParser(option_list = option_list)
    opt <- parse_args(opt_parser)
    return(opt)
}

combine_zns <- function(dirname) {
    fnames <- dir_ls(dirname, regexp = '.*[0-9]{6}.*')
    d_all <- map_dfr(fnames, read_delim, delim = ' ', 
                     col_types = cols()) %>%
             group_by(start) %>%
             filter(site_count == max(site_count)) %>%
             ungroup()
    return(d_all)
}

# parse args
opt <- args()
zns_all <- combine_zns(opt$directory)
write_delim(zns_all, path = opt$outfile, delim = ' ')
message('Complete.')
message(paste('File written to', opt$outfile))
