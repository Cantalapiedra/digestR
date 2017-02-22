#!/usr/bin/env Rscript
## fragment_size_selection
## 2017 CPCantalapiedra

library(data.table)

# Read command line arguments

args = commandArgs(trailingOnly=TRUE)

if (length(args)==3) {
  frag_file = args[1] # bed file
  min_size = strtoi(args[2])
  max_size = strtoi(args[3])
} else {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

# Read bed file

write(paste("Parsing bed file:", frag_file), stderr())

df_bed = fread(frag_file, header=TRUE, sep = "\t", verbose=FALSE, showProgress=FALSE)

write(paste("Rows:", nrow(df_bed)), stderr())

### Now perform size selection

write(paste("Size selection between", min_size, "and", max_size), stderr())

mins <- (df_bed$end - df_bed$start + 1)>=min_size
maxs <- (df_bed$end - df_bed$start + 1)<=max_size
filtered <- df_bed[mins & maxs,]

write(paste("Final rows:", nrow(filtered)), stderr())

# output

write.table(filtered, file = stdout(), row.names=FALSE, quote=FALSE, sep="\t")

## END