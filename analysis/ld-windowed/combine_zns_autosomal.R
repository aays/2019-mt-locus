# combine_zns_autosomal.R - combine zns files for a given chromosome
#
# usage:
# Rscript combine_zns_autosomal.R --directory [directory] --chrom [chrom] --outfile [file]

library(readr)
library(dplyr, warn.conflicts = FALSE)
library(purrr)
library(magrittr, warn.conflicts = FALSE)
library(fs)
library(optparse)

args <- function() {
    option_list <- list(
        make_option(c("-d", "--directory"), type = "character", default = NULL, 
                  help = "Directory containing ZnS files", metavar = "character"),
        make_option(c("-c", "--chrom"), type = "character", default = NULL, 
                  help = "Chromosome name", metavar = "character"),
        make_option(c("-o", "--outfile"), type = "character", default = NULL, 
                  help= "File to write to", metavar = "character")
    ) 
    opt_parser <- OptionParser(option_list = option_list)
    opt <- parse_args(opt_parser)
    return(opt)
}

combine_zns <- function(dirname, chrom) {
    fnames <- dir_ls(dirname, regexp = paste0(chrom, '[_\\:][0-9\\-]*\\.zns'))
    d_all <- map_dfr(fnames, read_delim, delim = ' ', 
                     col_types = cols())
    return(d_all)
}

# parse args
opt <- args()
zns_all <- combine_zns(opt$directory, opt$chrom)
write_delim(zns_all, path = opt$outfile, delim = ' ')
message('Complete.')
message(paste('File written to', opt$outfile))
