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

df_enz = read.table(enz_file, header=FALSE, sep = "\t", col.names = c("Enzyme", "Target"))

write.table(df_enz, stdout())

# Read genome from file

genome <- readDNAStringSet(genome_file, "fasta")

head(genome)

# Identify the MspI recongnition sites for each chromosomal entry
# Generate a dataframe with the length of MspI digested fragments
mdf=data.frame();
for (i in seq_along(Mmusculus)){
  #print(paste("Processing ",seqnames(Mmusculus)[i], sep=""))
  for (j in seq_along(Enzymes)){
    enz_target = Enzymes[i]$target
    m<-matchPattern(enz_target, Mmusculus[[i]])
    
    enz_target_len = seq_len(enz_target)
    starts<-start(gaps(m))-enz_target_len
    ends<-end(gaps(m))
    
    enz_name = Enzymes[i]$name
    temp_df<-data.frame(start=starts,end=ends,chr=seqnames(Mmusculus)[i],enzyme=enz_name) #actually end = ends
    temp_df$start<-replace(temp_df$start, temp_df$start == -3, 0)
    temp_df<-temp_df[c("chr","start","end","enzyme")]
    mdf<-rbind(mdf,temp_df)
}

# Extract the digested fragment length
# Keep only the fragments in the range of 40-600 bp
mdf$width=mdf$end-mdf$start
ml<-mdf[mdf$width>39&mdf$width<601,]
counts<-ddply(ml,.(width), nrow)
# Create a plot of the frequency of the fragment lengths (y-axis is logarithmic)
p<-ggplot(counts,aes(x=width, y=V1))+geom_line()
p+scale_y_continuous(trans=log2_trans())