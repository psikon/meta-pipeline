#!/usr/bin/env Rscript

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
message("Loading libary: blastr")
suppressPackageStartupMessages(library('blastr'))
message("Loading libary: ncbi")
suppressPackageStartupMessages(library('ncbi'))
message("Loading libary: metaR")
suppressPackageStartupMessages(library('metaR'))
suppressPackageStartupMessages(library('rmisc'))

message("Generate Taxonomy Database ...")
generate.TaxonomyReport(blast_db_path = opt$input,
                        metadata = list(SampleId=1),
                        coverage_threshold = opt$coverage,
                        min_match = 50,
                        taxon_db_path = opt$output,
                        bitscore_tolerance = opt$bitscore)
message("Finished")
