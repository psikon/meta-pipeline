#!/usr/bin/env Rscript

# load command line interface
suppressPackageStartupMessages(library('optparse'))

# specify the commandline interface
option_list <- list(
  make_option("--output", dest = "output", type = "character", 
              help = "path for taxonomy db"),
  make_option("--input", dest = "input", type ="character",
              help = "location of blastdb"),
  make_option("--bitscore_tolerance", dest = "bitscore", type = "numeric",
              default = 0.90, help = "bitscore_tolerance"),
  make_option("--coverage_threshold", dest = "coverage", type = "numeric",
              default = 0.30, help = "coverage_threshold"))
  make_option("--lib", dest = "lib", default = NULL, type = "character",
              help = "specify R libary position")

# init the commandline interface
opt <- parse_args(OptionParser(usage = "usage: %prog [options]",
                               option_list = option_list,
                               add_help_option = TRUE,
                               ))
# load needed libraries
message("Loading libary: blastr")
suppressPackageStartupMessages(library('blastr', lib.loc = opt$lib))
message("Loading libary: ncbi")
suppressPackageStartupMessages(library('ncbi', lib.loc = opt$lib))
message("Loading libary: metaR")
suppressPackageStartupMessages(library('metaR', lib.loc = opt$lib))
suppressPackageStartupMessages(library('rmisc', lib.loc = opt$lib))

message("Generate Taxonomy Database ...")
generate.TaxonomyReport(blast_db_path = opt$input,
                        metadata = list(SampleId=1),
                        coverage_threshold = opt$coverage,
                        min_match = 50,
                        taxon_db_path = opt$output,
                        bitscore_tolerance = opt$bitscore)
message("Finished")
