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
  seq <- genome[[i]]
  seq_name <- seq_names[[i]]
  #write(paste("Genome index:", i), stderr())
  
  ranges <- apply(df_enz, 1, f_enz, seq, seq_name)
  ranges <- do.call(rbind, ranges)
  return(ranges)
  #return()
}


write("Check sequences in genome...", stderr())
num_seqs <- length(genome)
write(paste("Num. sequences read:", num_seqs), stderr())

write(paste("Num. cores to balance tasks:", cores), stderr())

seqs_per_core = ceiling(num_seqs / cores)

write(paste("Num. sequences per core: ", seqs_per_core), stderr())
      
f_balance <- function(num_set, num_seqs, seqs_per_core, genome){
  #write(paste("Num set:", num_set), stderr())
  ini = 1+seqs_per_core*(num_set-1)
  end = ini+(seqs_per_core-1)
  if (end>num_seqs) end = num_seqs
  if (ini>num_seqs) return(NULL)
  
  #write(paste("Ini:", ini, "End:", end), stderr())
  return(genome[ini:end])
}

write("Balancing tasks...", stderr())
genome_list <- lapply(1:cores, f_balance, num_seqs, seqs_per_core, genome)
genome_list <- Filter(Negate(is.null), genome_list)
#genome_list

cores <- length(genome_list)
write(paste("Cores which will be used after balancing: ", cores), stderr())

# Function to find matches for each enzyme
# In theory, this could be done with matchPDict also avoiding
# the explicit loop (apply function)
f_enz <- function(enz, genome){
  
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

f_seq <- function(genome){
  #### Main loop: for each enzyme, look for targets  
  ranges <- apply(df_enz, 1, f_enz, genome)
  ranges <- do.call(rbind, ranges)
  return(ranges)
}

# Create cluster for parallel

write("Creating cluster...", stderr())
#cl <- makeCluster(cores, type="FORK", outfile = "debug.txt")
cl <- makeForkCluster(cores, outfile = "debug.txt")

write("\nLoop over genome list...", stderr())
loop_ranges <- parLapply(cl, genome_list, f_seq)
write("\nMerge results...", stderr())
all.ranges <- do.call(rbind, loop_ranges)

write("Releasing cluster...", stderr())
stopCluster(cl)

# Convert to data.table to sort the output date
write("Sorting results...", stderr())
all.ranges = data.table(all.ranges)
setorder(all.ranges, names, start, end, enzyme)

write("writing results...", stderr())

# Output data
write.table(all.ranges, file = stdout(), row.names=FALSE, quote=FALSE, sep="\t")

write("finished.", stderr())

q()

seq_names = names(genome)


#cl
#is.null(cl)

#write("Cluster export...", stderr())
#clusterExport(cl, "df_enz")
#clusterExport(cl, "f_enz")
#clusterExport(cl, "genome")





#q()



## END