#!/usr/bin/env Rscript
## fragment_type_selection
## 2017 CPCantalapiedra

library(data.table)

# Read command line arguments

args = commandArgs(trailingOnly=TRUE)

if (length(args)==2) {
  frag_file = args[1] # bed file
  frag_types_file = args[2]
} else {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

# Read fragment types file

df_frag_types <- fread(frag_types_file, header=FALSE, verbose=FALSE, showProgress=FALSE, sep="\t")

frag_types <- unlist(df_frag_types[,"V1"])
frag_types

# Read bed file

write(paste("Parsing bed file:", frag_file), stderr())

df_bed = fread(frag_file, header=TRUE, sep = "\t", verbose=FALSE, showProgress=FALSE)

head_bed <- head(df_bed, file=stderr())
write.table(head_bed, file = stderr(), row.names = FALSE, quote = FALSE, sep = "\t")
write(paste("Rows:", nrow(df_bed)), stderr())

### Filter by frag type

true_frag_types <- df_bed$fragtype %in% frag_types

filtered <- df_bed[true_frag_types,]#mat_lengths[mat_lengths$len>=min_size && mat_lengths$len>=max_size,]

write(paste("Final rows:", nrow(filtered)), stderr())

# output

write.table(filtered, file = stdout(), row.names=FALSE, quote=FALSE, sep="\t")

## END