#!/usr/bin/env Rscript
## digest_barplot
## 2017 CPCantalapiedra

library(data.table)

# Read command line arguments

args = commandArgs(trailingOnly=TRUE)

if (length(args)==2) {
  bed_file = args[1] # bed file
  out_file = args[2]
} else {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

# Read bed file

write(paste("Parsing bed file:", bed_file), stderr())

df_bed = fread(bed_file, header=TRUE, sep = "\t", verbose=FALSE, showProgress=FALSE)

#head_bed <- head(df_bed, file=stderr())
#write.table(head_bed, file = stderr(), row.names = FALSE, quote = FALSE, sep = "\t")

counts <- table(df_bed$names, df_bed$enzyme)

## Add totals to rows and columns

total.rows <- apply(counts, 1, function(x) sum(x))

counts <- cbind(counts, total=total.rows)

total.cols <- apply(counts, 2, function(x) sum(x))

counts <- rbind(counts, total=total.cols)

counts <- data.frame(counts)

# Sort by total
setorder(counts, -total)

# Transpose to order by column subtotals
t_counts <- data.frame(t(counts))

setorder(t_counts, -total)

# Transpose again to recover original table
counts <- data.frame(t(t_counts))

# Remove totals

counts <- counts[,!names(counts) %in% c("total")]
counts <- counts[rownames(counts)!="total",]

write("finished.", stderr())
counts

pdf(out_file)
barplot(as.matrix(counts))
dev.off()

## END