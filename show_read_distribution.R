#!/usr/bin/env Rscript

# load command line interface
suppressPackageStartupMessages(library('optparse'))

# specify the commandline interface
option_list <- list(
  make_option(c("-i", "--input"), dest = "datadir", type = "character",
              help = "output folder of meta-pipeline"),
  make_option(c("-r","--raw"), dest = "rawdir", type = "character",
              help = "folder with raw sequences"),
  make_option(c("-n", "--names"), dest = "samplenames", type = "character",
              help = "names of samples for x axis (in format '1,2,3,4,...')"),
  make_option(c("-o", "--output"), dest = "output", type = "character",
              default = "read_distribution", 
              help = "path and name for output"))

# init the commandline interface
opt <- parse_args(OptionParser(usage = "usage: %prog [options]",
                               option_list = option_list,
                               add_help_option = TRUE,
                               ))
# load needed libraries
suppressPackageStartupMessages(library('reshape2'))
suppressPackageStartupMessages(library('ggplot2'))
suppressPackageStartupMessages(library('scales'))
suppressPackageStartupMessages(library('ShortRead'))

# parse input from files
# find all files from Trimmomatic step
filtered_fastq <- list.files(opt$datadir, pattern = ".filtered.fastq", 
                             recursive = T, full.names = T)
# count sequences of filtered single end data
single_filtered <-  countLines(filtered_fastq[grep("*.single.filtered.fastq", filtered_fastq)])/4
# count sequences of filtered paired end data
filtered_fastq <- unique(countLines(filtered_fastq[-grep("*.single.filtered.fastq", filtered_fastq)]))/4
# count sequences of successful concatinated data
concat <- as.vector(countLines(list.files(opt$datadir, pattern = '*.extendedFrags.fastq', 
                                recursive = T, full.names = T))/4)
# count sequences of non duplicated data
no_dup <- as.vector(countLines(list.files(opt$datadir, pattern = '*.nodup.fasta', 
                                recursive = T, full.names = T))/2)
if(!is.null(opt$rawdir)) {
  # count sequences of raw data if sepcified
  raw = unique(countLines(list.files(opt$rawdir, pattern = '*.fastq', 
                              recursive = T, full.names = T))/4)
  data <- data.frame(raw, filtered_fastq, single_filtered,concat, no_dup)
} else {
  data <- data.frame(filtered_fastq,single_filtered,concat,no_dup)
}

# generate vector of rownames
samples <- unlist(strsplit(opt$samplenames,","))
if (length(samples)!= nrow(data)){
  # exit if new rownames unequal to number of rows
  message("Check samplenames: must be equal to number of samples")
  q("no", 1, FALSE)
} else {
  #adjust the rownames
  rownames(data) <- samples
}

# rearrange vector for ggplot2
data2 <- melt(as.data.frame(t(data)))
# new factor for legend 
if(!is.null(opt$rawdir)) {
  desc <- factor(rep(c("RAW", "paired end trimmed", "single end trimmed",
                       "concatinated", "single without duplicates"),
                     length=nrow(data2)),
                 levels=c("RAW", "paired end trimmed", "single end trimmed",
                          "concatinated", "single without duplicates"))
} else {
  desc <- factor(rep(c("paired end trimmed", "single end trimmed",
                       "concatinated", "single without duplicates"),
                     length=nrow(data2)),
                levels=c("paired end trimmed", "single end trimmed",
                         "concatinated", "single without duplicates"))
}
# combine to new data.frame
df <- cbind(data2, desc)
# generate output name
output <- paste0(opt$output,".pdf")
# build plot
pdf(output)
  ggplot(df, aes(x = variable, y = value, fill = desc)) + 
    # type of plot
    geom_bar(stat = "identity", position = "dodge") +
    # names of the axis
    xlab("\nSample") + ylab("Number of Reads")  +
    theme_bw() +
    # adjust Y-Axis
    scale_y_continuous(labels = comma, breaks = pretty_breaks(n=15)) +
    # change colors
    scale_fill_manual(values = c("brown2", "chartreuse4", "chartreuse1", 
                                "darkorchid1", "#6495ED")) + 
    # change legend title
    guides(fill = guide_legend("Type of reads")) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle("Distribution of Reads after preprocessing steps\n")
dev.off()
