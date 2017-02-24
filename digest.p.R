#!/usr/bin/env Rscript
## digestR
## 2017 CPCantalapiedra
## based on https://www.r-bloggers.com/restriction-digestion-of-eukaryotic-genomes-in-r/

library(Biostrings)
library(data.table)
library(parallel)

# Read command line arguments

args = commandArgs(trailingOnly=TRUE)

if (length(args)==3) {
  enz_file = args[1] # file with enzymes configuration
  genome_file = args[2] # file with input sequence used as target of those enzymes
  cores = strtoi(args[3])
} else if (length(args)==2) {
  enz_file = args[1] # file with enzymes configuration
  genome_file = args[2] # file with input sequence used as target of those enzymes
  cores = detectCores()
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
f_enz <- function(enz, subject, seq_name){
  #write("Current:", stderr())
  #write(enz, stderr())
  #write("Subject:", stderr())
  #write(seq_name, stderr())
  
  # search the RE pattern in the input sequence
  # fixed="subject" to avoid Ns while allowing UIPAC codes in enzyme target
  ranges <- matchPattern(enz["target"], subject, fixed="subject")
  
  num_matches = length(ranges)
  write(paste("\t", enz["name"], " - ", enz["target"], " --> Matches found: ", num_matches, sep = ""), file = stderr())
  
  # If there are hits, add the enzyme name and reorder columns
  if (num_matches > 0) {
    ranges <- data.frame(names=seq_name, start=start(ranges), end=end(ranges), enzyme=enz["name"], row.names = NULL)
  # when there are not hits just create an empty record
  } else {
    ranges <- data.frame(names=numeric(0), start=numeric(0), end=numeric(0), enzyme=numeric(0))
  }
  
  #write.table(ranges, stderr())
  
  return(ranges)
}

#### Main loop: for each enzyme, look for targets
#write("\nstart cutting...", stderr())
#all.ranges <- do.call(rbind, apply(df_enz, 1, f_enz, genome))

f_seq <- function(i, seq_names){
  #seq <- genome[[i]]
  #seq_name <- seq_names[[i]]
  #write(paste("Genome index:", i), stderr())
  
  #ranges <- apply(df_enz, 1, f_enz, seq, seq_name)
  #ranges <- do.call(rbind, ranges)
  #return(ranges)
  return()
}

write("Check sequences in genome...", stderr())
num_seqs <- length(genome)
seq_names = names(genome)

# Create cluster for parallel

write("Creating cluster...", stderr())
cl <- makeCluster(cores, type="FORK", outfile = "debug.txt")

write("Cluster export...", stderr())
clusterExport(cl, "df_enz")
clusterExport(cl, "f_enz")
clusterExport(cl, "genome")

write("\nLoop over genome...", stderr())
loop_ranges <- parLapply(cl, 1:num_seqs, f_seq, seq_names)
all.ranges <- do.call(rbind, loop_ranges)

stopCluster(cl)
#q()

# Convert to data.table to sort the output date
all.ranges = data.table(all.ranges)
setorder(all.ranges, names, start, end, enzyme)

write("writing results...", stderr())

# Output data
write.table(all.ranges, file = stdout(), row.names=FALSE, quote=FALSE, sep="\t")

write("finished.", stderr())

## END