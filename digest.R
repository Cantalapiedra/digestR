#!/usr/bin/env Rscript
## digest.R
## 2017 CPCantalapiedra
## based on https://www.r-bloggers.com/restriction-digestion-of-eukaryotic-genomes-in-r/

library(Biostrings)
library(data.table)

# Read command line arguments

args = commandArgs(trailingOnly=TRUE)

if (length(args)==2) {
  enz_file = args[1] # file with enzymes configuration
  genome_file = args[2] # file with input sequence used as target of those enzymes
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
#write(paste("\t", names(genome), "of length:", width(genome)), stderr())

# Function to find matches for each enzyme
# In theory, this could be done with matchPDict also avoiding
# the explicit loop (apply function)
f_enz <- function(enz){
  
  # search the RE pattern in the input sequence
  # fixed="subject" to avoid Ns while allowing UIPAC codes in enzyme target
  match <- vmatchPattern(enz["target"], genome, fixed="subject")
  ranges <- unlist(match)
  
  write(paste("\t", enz["name"], " - ", enz["target"], " --> Matches found: ", length(ranges), sep = ""), file = stderr())
  
  # If there are hits, add the enzyme name and reorder columns
  if (length(ranges)>0) {
    ranges <- as.data.frame(ranges)
    ranges <- cbind(ranges, enzyme=enz["name"], row.names = NULL)
    ranges <- ranges[,c("names", "start", "end", "enzyme")]
    #write.table(ranges, stderr())
  } else { # when there are not hits just create an empty record
    ranges <- as.data.frame(ranges)
    ranges <- cbind(ranges, enzyme=numeric(0), row.names = NULL)
    ranges <- ranges[,c("names", "start", "end", "enzyme")]
  }
  
  return(ranges)
}

write("\nstart cutting...", stderr())

#### Main loop: for each enzyme, look for targets
all.ranges <- do.call(rbind, apply(df_enz, 1, f_enz))

# Convert to data.table to sort the output date
all.ranges = data.table(all.ranges)
setorder(all.ranges, names, start, end, enzyme)

write("writing results...", stderr())

# Output data
write.table(all.ranges, file = stdout(), row.names=FALSE, quote=FALSE, sep="\t")

write("finished.", stderr())

## END