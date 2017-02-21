#!/usr/bin/env Rscript
## digestR
## 2017 CPCantalapiedra
## based on https://www.r-bloggers.com/restriction-digestion-of-eukaryotic-genomes-in-r/

### Load the needed libraries
# Biostrings is needed for pattern identification
# plyr and reshape2 are needed for manipulating the data format

library(Biostrings)
#library(BSgenome.Mmusculus.UCSC.mm10) # HOW TO LOAD THE Hvulg genome
library(plyr)
library(ggplot2)
library(reshape2)
library(scales)

# Read command line arguments

args = commandArgs(trailingOnly=TRUE)

if (length(args)==2) {
  # default output file
  enz_file = args[1]
  genome_file = args[2]
  
  write(enz_file, stdout())
  write(genome_file, stdout())
  
} else {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

# Read enzymes from a configuration file

df_enz = read.table(enz_file, header=FALSE, sep = "\t", col.names = c("name", "target"))

write.table(df_enz, stdout())

# Read genome from file

genome <- readDNAStringSet(genome_file, "fasta")

genome

# Identify the MspI recongnition sites for each chromosomal entry
# Generate a dataframe with the length of MspI digested fragments
mdf=data.frame();

f_enz <- function(enz){
  write(paste("\t", enz["name"], "\t", enz["target"]), stdout())
  #count <- vcountPattern(enz["target"], genome)
  #write(count, stdout())
  match <- vmatchPattern(enz["target"], genome, fix=FALSE)
  write.table(match, stdout())
}

write("START", stdout())
apply(df_enz, 1, f_enz)

write("END", stdout())

## END