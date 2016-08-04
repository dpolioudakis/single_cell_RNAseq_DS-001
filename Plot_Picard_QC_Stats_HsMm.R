# Damon Polioudakis
# 2016-06-30
# Process Picard QC statistics from RNAseq

# Fragments mapped by location must be divided by correct read size
################################################################################

rm(list=ls())
sessionInfo()
library(ggplot2)
library(reshape2)

# Load data and assign variables

# alignment_stats.txt
picAlign <- read.table("../data/QC/human_mouse/SxaQSEQsXap108L1/alignment_stats.txt"
                       , row.names = 1, skip = 1)
# remove unused column names from alignment_stats.txt
cNames <- unlist(strsplit(readLines(con = "../data/QC/human_mouse/SxaQSEQsXap108L1/alignment_stats.txt"
                             , n = 1), split = " "))
colnames(picAlign) <- cNames[ -c(1, (length(cNames)-2):length(cNames))]

# rnaseq_stats.txt
picSeq <- read.table("../data/QC/human_mouse/SxaQSEQsXap108L1/rnaseq_stats.txt", skip = 1, row.names = 1)
# remove unused column names from alignment_stats.txt
# (Picard outputs column names for statistics even if they are not calculated)
cNames <- unlist(strsplit(readLines(con = "../data/QC/human_mouse/SxaQSEQsXap108L1/rnaseq_stats.txt"
                                    , n = 1), split = " "))
colnames(picSeq) <- cNames[!cNames %in% c("RIBOSOMAL_BASES"
                  , "PCT_RIBOSOMAL_BASES", "SAMPLE", "LIBRARY", "READ_GROUP")]

# rnaseq_stats_Transcript_Coverage.txt
# Add NA to blanks - I think these are from samples with almost no coverage
txCovDF <- read.table("../data/QC/human_mouse/SxaQSEQsXap108L1/rnaseq_stats_Transcript_Coverage.txt"
                      , fill = TRUE)

# gcbias_summary.txt
gcBiasSummaryDF <- read.table("../data/QC/human_mouse/SxaQSEQsXap108L1/gcbias_summary.txt", header = TRUE)

# gcbias_stats.txt
gcBiasStatsDF <- read.table("../data/QC/human_mouse/SxaQSEQsXap108L1/gcbias_stats.txt", header = FALSE)

# duplication_stats.txt
# Add NA to blanks
dupDF <- read.table("../data/QC/human_mouse/SxaQSEQsXap108L1/duplication_stats.txt"
                    , fill = TRUE, header = TRUE)
# Shift columns so they align with headers
dupDF[ ,c(3:ncol(dupDF))] <- dupDF[ ,c(4:ncol(dupDF))]

# Make out graphs directory
dir.create("../analysis/graphs", recursive = TRUE)

# Set ggplot2 theme
theme_set(theme_bw())
theme_set(theme_get() + theme(text=element_text(size=18)))
theme_update(plot.title = element_text(size = 16))
################################################################################

# RNAseq Percentage of Bases Aligned by Location

pctHq <- picAlign[,"PF_HQ_ALIGNED_BASES"]/picAlign[,"PF_ALIGNED_BASES"]

ggDF <- cbind(pctHq, picSeq[,c(14, 10:13)])
colnames(ggDF) <- c("High\nQuality"
                    , "mRNA\n"
                    , "Protein\nCoding"
                    , "Untranslated\nRegion"
                    , "Intronic\nRegion"
                    , "Intergenic\nRegion")
ggDF <- melt(ggDF)

ggplot(ggDF, aes(x = variable, y = value, col = variable)) +
  geom_boxplot() +
  guides(col = FALSE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Percent of Aligned Bases") +
  xlab("Location") +
  ggtitle(paste0("Plot_Picard_QC_Stats_HsMm.R"
                 , "\nRNA-seq QC metrics across samples"
                 , "\nPercent of Aligned Bases by Location"
                 , "\n"))
ggsave("../analysis/graphs/Plot_Picard_QC_HsMm_Percent_Location.pdf"
       , width = 7, height = 7)
################################################################################

# Boxplot number of fragments mapping to different types of location
# Change reads to fragments for paired end data

# Divide by read size
ggDF <- picSeq[, c("PF_BASES", "PF_ALIGNED_BASES", "CODING_BASES", "UTR_BASES", "INTRONIC_BASES", "INTERGENIC_BASES")]/64
ggDF <- cbind(ggDF$PF_BASES, ggDF$PF_ALIGNED_BASES, mRNA = (ggDF$CODING_BASES + ggDF$UTR_BASES)
              , ggDF[, c("CODING_BASES", "UTR_BASES", "INTRONIC_BASES", "INTERGENIC_BASES")])
colnames(ggDF) <- c("Total Reads", "Aligned", "mRNA", "Protein Coding", "UTR", "Intronic", "Intergenic")
ggDF <- melt(ggDF)

ggplot(ggDF, aes(x = variable, y = value, col = variable)) +
  geom_boxplot(outlier.shape = NA) +
  # Adjust limits after outlier removal
  coord_cartesian(ylim = range(boxplot(ggDF[ggDF$variable == "Total Reads", ]$value
                                       , plot = FALSE)$stats) * c(0, 1.1)) +
  guides(col = FALSE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Number of Mapped Reads") +
  xlab("Location") +
  ggtitle(paste0("Plot_Picard_QC_Stats_HsMm.R"
                 , "\nNumber of Mapped Reads by Location"
                 , "\n"))
ggsave("../analysis/graphs/Plot_Picard_QC_HsMm_Aligned_Location.pdf"
       , width = 7, height = 7)
################################################################################

# Transcript Coverage

keep <- match(unique(txCovDF[ ,1]), txCovDF[ ,1])
txCovDF <- txCovDF[keep, ]
rownames(txCovDF) <- txCovDF[ ,1]
txCovDF.cols <- as.character(txCovDF[1, seq(4, 204, by = 2)])
txCovDF <- txCovDF[ ,seq(5, 205, by = 2)]
colnames(txCovDF) <- txCovDF.cols
txCovqtlDF <- apply(txCovDF, 2, quantile, c(0.025, 0.5, 0.975), na.rm = TRUE)

pdf("../analysis/graphs/Plot_Picard_QC_HsMm_Transcript_Coverage.pdf"
    , width = 7, height = 7)
plot(x = as.numeric(colnames(txCovqtlDF)), y = txCovqtlDF[2,]
     , xlab = "Percentile of gene body (5' -> 3')"
     , ylab = "Coverage Relative to Whole Transcript"
     , pch = 19, ylim = c(0, 7), cex.lab = 1.4
     , main = "Plot_Picard_QC_Stats_HsMm.R
Relative transcript coverage - median with 95% CIs across samples")
lines(x = as.numeric(colnames(txCovqtlDF)), y = txCovqtlDF[2,], col = "black")
lines(x = as.numeric(colnames(txCovqtlDF)), y = txCovqtlDF[1,], col = "grey")
lines(x = as.numeric(colnames(txCovqtlDF)), y = txCovqtlDF[3,], col = "grey")
dev.off()
################################################################################

## GC Bias

# Normalized coverage by GC %
gcBias <- matrix(NA, nrow = nrow(gcBiasStatsDF), ncol = 101)
colnames(gcBias) <- as.character(gcBiasStatsDF[1, seq(2, 602, by = 6)])
rownames(gcBias) <- gcBiasStatsDF[ ,1]
gcBias <- gcBiasStatsDF[ ,seq(6, 607, by = 6)]
gc.coverage.quant <- apply(gcBias, 2, quantile, c(0.025, 0.5, 0.975))
colnames(gc.coverage.quant) <- as.character(gcBiasStatsDF[1, seq(2, 602, by = 6)])

# Read quality score by GC %
gcBiasqual <- matrix(NA, nrow = nrow(gcBiasStatsDF), ncol = 101)
colnames(gcBiasqual) <- gcBiasStatsDF[1, seq(2, 602, by = 6)]
rownames(gcBiasqual) <- gcBiasStatsDF[ ,1]
gcBiasqual <- gcBiasStatsDF[ ,seq(5, 607, by = 6)]
gc.qual.quant <- apply(gcBiasqual, 2, quantile, c(0.025, 0.5, 0.975))

# Proportion of 100bp bins
gcBiasbin <- matrix(NA, nrow = nrow(gcBiasStatsDF), ncol = 101)
colnames(gcBiasbin) <- gcBiasStatsDF[1, seq(2, 602, by = 6)]
rownames(gcBiasbin) <- gcBiasStatsDF[ ,1]
gcBiasbin <- gcBiasStatsDF[ ,seq(3, 607, by = 6)]
gc.bins.quant <- apply(gcBiasbin, 2, quantile, c(0.5))
gc.bins.quant <- gc.bins.quant / sum(gc.bins.quant)

# Plots
pdf("../analysis/graphs/Plot_Picard_QC_HsMm_GC_Bias_Stats.pdf"
    , width = 7, height = 7)
# Normalized coverage by GC %
par(mfrow = c(3,1))
plot(x = as.numeric(colnames(gc.coverage.quant))
     , y = gc.coverage.quant[2, ]
     , xlab = "% GC content"
     , ylab = "Normalized coverage"
     , pch = 19, cex.lab = 1.4, cex.axis = 1.2
     , ylim = c(0, 30)
     , main = "Coverage by GC percentage - median with 95% CIs across samples")
lines(x = as.numeric(colnames(gc.coverage.quant))
      , y = gc.coverage.quant[2, ], col = "black")
lines(x = as.numeric(colnames(gc.coverage.quant))
      , y = gc.coverage.quant[1, ], col = "grey")
lines(x = as.numeric(colnames(gc.coverage.quant))
      , y = gc.coverage.quant[3,], col = "grey")

# Read quality score by GC %
plot(x = as.numeric(colnames(gc.coverage.quant))
     , y = gc.qual.quant[2, ]
     , xlab = "%GC content"
     , ylab = "Read quality score"
     , pch = 19, cex.lab = 1.4, cex.axis = 1.2
     , ylim = c(0, 36)
     , main = "Read quality by GC percent - median with 95% CIs across samples")
lines(x = as.numeric(colnames(gc.coverage.quant))
      , y = gc.qual.quant[2,], col = "black")
lines(x = as.numeric(colnames(gc.coverage.quant))
      , y = gc.qual.quant[1,], col = "grey")
lines(x = as.numeric(colnames(gc.coverage.quant))
      , y = gc.qual.quant[3,], col = "grey")

plot(density(x = as.numeric(colnames(gc.coverage.quant)), y = gc.bins.quant)
     , xlab = "%GC content"
     , ylab = "Proportion of 100bp bins"
     , pch = 19, cex.lab = 1.4, cex.axis = 1.2
     , ylim = c(0, max(gc.bins.quant) + 0.01)
     , main = "Proportion of bins corresponding to each GC percentile")
dev.off()
################################################################################

### Duplication Metrics

## Percent Duplicate for all capture sites

# For Single End RNAseq
ggDF <- data.frame(PERCENT_DUPLICATION = dupDF$PERCENT_DUPLICATION * 100
                   , SAMPLE = gsub("/Picard_Stats", "", dupDF$SAMPLE))
# Order by percent duplication
# ggDF <- ggDF[order(ggDF$PERCENT_DUPLICATION), ]
ggDF$SAMPLE <- factor(ggDF$SAMPLE, levels = ggDF$SAMPLE)
mPctDup <- signif(mean(ggDF$PERCENT_DUPLICATION), 4)

ggplot(ggDF, aes(x = SAMPLE, y = PERCENT_DUPLICATION)) +
  geom_bar(stat = "identity") +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  ylab("Percentage Duplication") +
  xlab("Samples") +
  ggtitle(paste0("Plot_Picard_QC_Stats_HsMm.R"
                , "\nPercentage of Mapped Sequence Marked as Duplicate"
                , "\nMean Percent Duplicate: ", mPctDup))
ggsave("../analysis/graphs/Plot_Picard_QC_HsMm_Percent_Duplication.pdf"
       , width = 7, height = 7)
################################################################################

# Compile QC metrics

## Compile most important QC metrics - choose the metrics that you think are valuable to the analysis and to sequencing statistics
## Metrics concering read depth (ie. Total Reads) and GC bias (ie. 5'/3' bias) should always be included

PQCdat <- cbind(picAlign[ ,c(2, 6, 9, 21)], picSeq, GC_DROPOUT = gcBiasSummaryDF[ ,c(6)]
                , dupDF[ ,c(7, 8)])
write.csv(PQCdat,"../metadata/PicardToolsQC.csv")
dev.off()
