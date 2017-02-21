#!/usr/bin/env Rscript
## digestR
## 2017 CPCantalapiedra
## based on https://www.r-bloggers.com/restriction-digestion-of-eukaryotic-genomes-in-r/

library(Biostrings)
library(data.table)

# Read command line arguments

args = commandArgs(trailingOnly=TRUE)

if (length(args)==2) {
  # default output file
  enz_file = args[1]
  genome_file = args[2]
  
} else {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

# Read enzymes from a configuration file

write(paste("Parsing enzymes configuration file:", enz_file), stderr())

df_enz = read.table(enz_file, header=FALSE, sep = "\t", col.names = c("name", "target"))

write.table(df_enz, file = stderr(), row.names = FALSE, quote = FALSE, sep = "\t")

# Read genome from file

write(paste("\nParsing fasta input file:", genome_file), stderr())
genome <- readDNAStringSet(genome_file, "fasta")

write(paste("Num. sequences read:", length(genome)), stderr())
write(paste("\t", names(genome), "of length:", width(genome)), stderr())

# Function to find matches for each enzyme
# In theory, this could be done with matchPDict also avoiding
# the explicit loop (apply function)
f_enz <- function(enz){
  
  match <- vmatchPattern(enz["target"], genome, fix=FALSE)
  ranges <- unlist(match)
  
  write(paste("\t", enz["name"], " - ", enz["target"], " --> Matches found: ", length(ranges), sep = ""), file = stderr())
  
  if (length(ranges)>0) {
    ranges <- as.data.frame(ranges)
    ranges <- cbind(ranges, enzyme=enz["name"], row.names = NULL)
    ranges <- ranges[,c("names", "start", "end", "enzyme")]
    #write.table(ranges, stderr())
  } else {
    ranges <- as.data.frame(ranges)
    ranges <- cbind(ranges, enzyme=numeric(0), row.names = NULL)
    ranges <- ranges[,c("names", "start", "end", "enzyme")]
  }
  
  return(ranges)
}

write("\nstart cutting...", stderr())

all.ranges <- do.call(rbind, apply(df_enz, 1, f_enz))
all.ranges = data.table(all.ranges)
setorder(all.ranges, names, start, end, enzyme)

write("writing results...", stderr())

write.table(all.ranges, file = stdout(), row.names=FALSE, quote=FALSE, sep="\t")

write("finished.", stderr())

## END