# Damon Polioudakis
# 2016-06-30
# Process Picard QC statistics from RNAseq

# Plot statistics from Drop-seq adaptor trimming, poly-A trimming, cell barcode
# filter, and molecular barcode filter
################################################################################

rm(list=ls())
sessionInfo()

# Load data and assign variables

dirs <- list.files("../data/bam/human_mouse/SxaQSEQsXap108L1")

# dirs <- c("N702", "N703")

dir.create("../analysis/graphs", recursive = TRUE)
################################################################################

# Adaptor trimming stats
pdf("../analysis/graphs/Plot_Dropseq_Processing_Stats_Adapter_Trimming.pdf")
for (dir in dirs) {
  df <- read.table(paste0("../data/bam/human_mouse/SxaQSEQsXap108L1/", dir
                          , "/adapter_trimming_report.txt")
                   , header = TRUE, stringsAsFactors=F)
  mn <- df[1,1]
  df <- data.frame(df[-c(1:2), ])
  colnames(df) <- c("BIN", "VALUE")
  df$VALUE <- as.numeric(df$VALUE)
  print(head(df))
  nTrimmed <- sum(df$VALUE)
  
  barplot(df$VALUE, col = "blue"
          , names.arg = df$BIN
          , main = paste0("Dropseq Adapter Trimming - Histogram: ", dir
                          ,"\nMean: ", mn
                          ,"\nNumber Trimmed: ", nTrimmed)
          , xlab = "SMART adapter bases trimmed"
          , ylab = "Frequency"
  )
}
dev.off()

# Poly-A trimming stats
pdf("../analysis/graphs/Plot_Dropseq_Processing_Stats_PolyA_Trimming.pdf")
for (dir in dirs) {
  df <- read.table(paste0("../data/bam/human_mouse/SxaQSEQsXap108L1/", dir
                          , "/polyA_trimming_report.txt")
                   , header = TRUE, stringsAsFactors=F)
  mn <- df[1,1]
  df <- data.frame(df[-c(1:2), ])
  colnames(df) <- c("BIN", "VALUE")
  df$VALUE <- as.numeric(df$VALUE)
  print(head(df))
  nTrimmed <- sum(df$VALUE)
  
  barplot(df$VALUE, col = "blue"
          , names.arg = df$BIN
          , main = paste0("Dropseq PolyA Trimming: ", dir
                          ,"\nMean: ", mn
                          ,"\nNumber Trimmed: ", nTrimmed)
          , xlab = "PolyA bases trimmed"
          , ylab = "Frequency"
  )
}
dev.off()

# Cell barcode filter stats
pdf("../analysis/graphs/Plot_Dropseq_Processing_Stats_CellBarcode.pdf")
for (dir in dirs) {
  df <- read.table(paste0("../data/bam/human_mouse/SxaQSEQsXap108L1/", dir
                          , "/unaligned_tagged_Cellular.bam_summary.txt")
                   , header = TRUE, stringsAsFactors=F)
  df$VALUE <- as.numeric(df$num_barcodes)
  print(head(df))
  pctFail <- round((sum(df$VALUE[-1])/df$VALUE[1])*100, 2)
  
  barplot(df$VALUE, col = "blue"
          , names.arg = df$num_failed_bases
          , main = paste0("Cellular Barcode Failed Bases: ", dir
          , "\nPercent reads filtered from failed cellular barcode: ", pctFail)
          , xlab = "Number of bases below map quality 10"
          , ylab = "Frequency"
  )
}
dev.off()

# Molecular barcode filter stats
pdf("../analysis/graphs/Plot_Dropseq_Processing_Stats_MolecularBarcode.pdf")
for (dir in dirs) {
  df <- read.table(paste0("../data/bam/human_mouse/SxaQSEQsXap108L1/", dir
                          , "/unaligned_tagged_Molecular.bam_summary.txt")
                   , header = TRUE, stringsAsFactors=F)
  df$VALUE <- as.numeric(df$num_barcodes)
  print(head(df))
  pctFail <- round((sum(df$VALUE[-1])/df$VALUE[1])*100, 2)
  
  barplot(df$VALUE, col = "blue"
          , names.arg = df$num_failed_bases
          , main = paste0("Molecular Barcode Failed Bases: ", dir
          , "\nPercent reads filtered from failed molecular barcode: ", pctFail)
          , xlab = "Number of bases below map quality 10"
          , ylab = "Frequency"
  )
}
dev.off()