#!/usr/bin/env Rscript
## fragment_len_dist
## 2017 CPCantalapiedra

library(data.table)

# Read command line arguments

args = commandArgs(trailingOnly=TRUE)

if (length(args)==2) {
  frag_file = args[1] # bed file
  out_file = args[2]
} else {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

# Read bed file

write(paste("Parsing bed file:", frag_file), stderr())

df_bed = fread(frag_file, header=TRUE, sep = "\t", verbose=FALSE, showProgress=FALSE)

head_bed <- head(df_bed, file=stderr())
write.table(head_bed, file = stderr(), row.names = FALSE, quote = FALSE, sep = "\t")
write(paste("Rows:", nrow(df_bed)), stderr())

# Subsample to speed up calculations and plotting

if (nrow(df_bed)>1000000){
  sample_size = floor(nrow(df_bed)*1/100)
  write(paste("Sampling", sample_size, "elements..."), stderr())
  df_bed <- df_bed[sample(nrow(df_bed), size = sample_size),]
}

# Obtain fragment types

frag_types <- unique(df_bed$"fragtype")
write("Fragment types:", stderr())
frag_types

### Now calculate length distributions

write("Calculating lengths...", stderr())
mat_lengths <- cbind(df_bed, len=(df_bed$end - df_bed$start + 1))
head(mat_lengths)

# Create histograms

write("Creating histograms...", stderr())

create_hist <- function(frag_type, mat_lengths){
  write(paste("Next fragment type: ", frag_type), stderr())
  
  frag_type_data = mat_lengths[mat_lengths$fragtype==frag_type,]
  lens = unlist(frag_type_data[,"len"])
  h = hist(lens, plot=FALSE, breaks=400)
  h$density = h$counts/sum(h$counts)*100
  plot(h, freq=FALSE, main=frag_type, xlim=c(0, 2000), xlab="fragment size", ylab="percentage")
}

pdf(out_file)

# histograms
#frag_types <- frag_types[frag_types %in% "EcoRI:EcoRI"]
dummy <- lapply(frag_types, create_hist, mat_lengths)

# boxplot
write("Creating boxplot...", stderr())

#mat_lengths <- mat_lengths[mat_lengths$fragtype %in% c("MseI:EcoRI", "EcoRI:EcoRI"),] 
boxplot(len~fragtype, data=mat_lengths, ylim=c(0, 2000))

invisible(dev.off())

#warnings()

## END