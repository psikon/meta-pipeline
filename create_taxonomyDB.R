!/usr/bin/env Rscript

# load command line interface
suppressPackageStartupMessages(library('optparse'))

# specify the commandline interface
option_list <- list(
  make_option(c("-o", "--output"), dest = "output", type = "character", 
              help = "path for taxonomy db"))
  make_option(c("-i","--input"), dest = "input", type ="character",
              help = "location of blastdb")
  make_option(c("--bitscore_tolerance", ), dest = "bitscore", type = "numeric",
              default = 0.90, help = "bitscore_tolerance"),
  make_option(c("--coverage_threshold"), dest = "coverage", type = "numeric",
              default = 0.30, help = "coverage_threshold"),

# init the commandline interface
opt <- parse_args(OptionParser(usage = "usage: %prog [options]",
                               option_list = option_list,
                               add_help_option = TRUE,
                               ))
# load needed libraries
suppressPackageStartupMessages(library('blastr'))
suppressPackageStartupMessages(library('ncbi'))
suppressPackageStartupMessages(library('metaR'))
suppressPackageStartupMessages(library('rmisc'))
