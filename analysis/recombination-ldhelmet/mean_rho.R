# mean_rho.R - returns mean rho across input ldhelmet files
#

library(readr)
library(dplyr, warn.conflicts = FALSE)
library(purrr)
library(magrittr)
library(optparse)

args <- function() {
    option_list <- list(
        make_option(c("-d", "--dir"), type = "character", default = NULL, 
                  help = "Directory containing LDhelmet files", metavar = "character"),
        make_option(c("-o", "--outfile"), type = "character", default = NULL, 
                  help= "Output file name", metavar = "character")
    ) 
    opt_parser <- OptionParser(option_list = option_list)
    opt <- parse_args(opt_parser)
    return(opt)
}

weighted_mean <- function(df) {
    output <- df %>%
        transmute(length = right_snp - left_snp,
                  weighted_rho = mean * length) %>%
        summarise_all(sum) %>%
        mutate(rho_per_bp = weighted_rho / length)
    return(output)
}

# read in files from dir

opt <- args()
ldhelmet_cols <- c('left_snp', 'right_snp', 'mean', 'p025', 'p50', 'p975')
ldhelmet_files <- map(list.files(opt$dir, full.names = TRUE),
                      read_delim, delim = ' ', skip = 3,
                      col_types = cols(), col_names = ldhelmet_cols)
output <- map_dfr(ldhelmet_files, weighted_mean, .id = 'name')

write_delim(output, path = opt$outfile, delim = ' ')

