# Damon Polioudakis
# 2016-08-01
# Plot Y chromosome and Xist expression for Drop-seq DS-001
################################################################################

rm(list=ls())
sessionInfo()

require(reshape2)
require(ggplot2)
require(biomaRt)

## Load data and assign variables

# dirs <- c("N702", "N703")
dirs <- list.files("../data/digital_gene_expression/Homo_sapiens.GRCh37.75.dna.primary_assembly_NoERCC/SxaQSEQsXap108L1/")
# dir <- "N702"

inParent <- "../data/digital_gene_expression"

dir.create("../analysis/graphs", recursive = TRUE)
dir.create("../analysis/tables/", recursive = TRUE)

tReadsDF <- read.table("../data/QC/Reads_Per_Sample.txt", fill = TRUE
                       , header = TRUE)

outGraph <- "../analysis/graphs/Ychr_Xist_Expression_"
graphTitle <- "Ychr_Xist_Expression.R"

# Set ggplot2 theme
theme_set(theme_bw())
theme_set(theme_get() + theme(text=element_text(size=18)))
theme_update(plot.title = element_text(size = 16))
################################################################################

### Functions

## Convert Ensembl IDs to Gene Symbols

AddChromosomeHuman <- function (ensemblList) {
  ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
  moduleGenes <- data.frame(ensemblList)
  # bioMart manual:
  #http://bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/biomaRt.pdf
  # Attribute: values to retrieve
  # Filters: input query type
  # Values: input query
  #ensembl <- useMart("ensembl")
  ensembl <- useMart("ENSEMBL_MART_ENSEMBL", host="www.ensembl.org")
  ensembl <- useDataset("hsapiens_gene_ensembl", mart=ensembl)
  # Data frame of module Ensembl IDs and gene symbols
  ensemblGeneSymDF <- getBM(  attributes = c("ensembl_gene_id", "hgnc_symbol", "chromosome_name")
    , filters = "hgnc_symbol"
    , values = moduleGenes
    , mart = ensembl
  )
  ensemblGeneSymDF
}
################################################################################

### Filter digital gene expression tables for > 10,000 human and < 1000 mouse reads

# Convert FASTQ line count to read count (divide by 4)
tReadsDF$TOTAL_FRAGMENTS <- tReadsDF$LINES_IN_FASTQ / 4

# Loop through sample directories, make list of filtered digital gene expression
# data frames
dgeLDF <- list()
for (dir in dirs) {
  
  # Load human counts pre UMI collapse from human mouse mapping to filter
  hsDat <- read.table(paste0(inParent, "/human_mouse/SxaQSEQsXap108L1/", dir
                             , "/out_gene_exon_tagged_human.counts.summary.txt"), header = TRUE)
  
  print("")
  print(dir)
  str(hsDat)
  # List of barcodes with > X human exonic reads
  print(dim(hsDat[hsDat$NUM_TRANSCRIPTS > 10000, ]))
  hsBarcodes <- hsDat[hsDat$NUM_TRANSCRIPTS > 10000, ]$CELL_BARCODE
  
  # Load mouse counts pre UMI collapse from human mouse mapping to filter
  mmDat <- read.table(paste0(inParent, "/human_mouse/SxaQSEQsXap108L1/", dir
                             , "/out_gene_exon_tagged_mouse.counts.summary.txt"), header = TRUE)
  str(mmDat)
  # List of barcodes with < X mouse exonic reads
  print(dim(mmDat[mmDat$NUM_TRANSCRIPTS < 1000, ]))
  mmBarcodes <- mmDat[mmDat$NUM_TRANSCRIPTS < 1000, ]$CELL_BARCODE
  
  hsDgeDF <- read.table(paste0(inParent
                                     , "/Homo_sapiens.GRCh37.75.dna.primary_assembly_NoERCC/SxaQSEQsXap108L1/"
                                     , dir, "/out_gene_exon_tagged_human.dge.txt.gz")
                              , header = TRUE, row.names = 1)
  
  # Filtering digital gene expression table to list of cell barcodes with > X
  # human and < X mouse exonic reads
  hsFtDgeDF <- hsDgeDF[ ,colnames(hsDgeDF) %in% mmBarcodes & colnames(hsDgeDF) %in% hsBarcodes]
  dim(hsFtDgeDF)
  dgeLDF[[dir]] <- hsFtDgeDF
}
################################################################################

### Plot Xist and Y chromosome expression

## Convert Ensembl IDs to Gene Symbols

for (sample in names(dgeLDF)) {
  dgeDF <- dgeLDF[[sample]]
  # Dataframe of ensembl IDs and chromosome
  # Human
  hsEnsemblChrDF <- AddChromosomeHuman(row.names(dgeDF))
  print(dim(hsEnsemblChrDF))
  
  # Add chromosome to expression dataframe
  # Human
  dgeLDF[[sample]] <- merge(hsEnsemblChrDF, dgeDF, by.x = "hgnc_symbol", by.y = "row.names")
}

## Y expression
pdf(paste0(outGraph, "Mean_Ychr_Expression.pdf"))
for (sample in names(dgeLDF)) {
  dgeDF <- dgeLDF[[sample]]
  yHsFtExDF <- dgeDF[dgeDF$chromosome_name == "Y", ]
  barplot(sort(apply(yHsFtExDF[ ,-c(1:3)], 2, mean))
          , xaxt = "n", xlab = "Cells"
          , ylab = "Expression (counts)"
          , main = paste0(graphTitle
                          , "\n", sample, ": Drop-seq mean Y chromosome expression"
                          , "\n", "> 10^4 human reads aligned, < 10^3 mouse reads aligned"))
}
dev.off()

## X vs Y expression
pdf(paste0(outGraph, "Ychr_Vs_Xchr_Expression.pdf"))
for (sample in names(dgeLDF)) {
  dgeDF <- dgeLDF[[sample]]
  xyDF <- dgeDF[dgeDF$chromosome_name == "Y" | dgeDF$chromosome_name == "X", ]
  mnDF <- aggregate(xyDF[ ,-c(1:3)], list(chromosome = xyDF$chromosome_name), sum)
  row.names(mnDF) <- mnDF$chromosome
  mnDF <- mnDF[ ,-1]
  ggDF <- data.frame(t(mnDF))
  ggDF <- ggDF[order(ggDF$Y), ]
  graph <- ggplot(ggDF, aes(x = X, y = Y)) +
    geom_point() +
    ylab("Y chromosome expression (total counts)") +
    xlab("X chromosome expression (total counts)") + 
    ggtitle(paste0(graphTitle
      , "\n", sample, ": Drop-seq X versus Y chromosome expression"
      , "\n", "> 10^4 human reads aligned, < 10^3 mouse reads aligned"
      , "\n"))
  print(graph)
}
dev.off()

## Xist expression
pdf(paste0(outGraph, "Xist_Expression.pdf"))
for (sample in names(dgeLDF)) {
  dgeDF <- dgeLDF[[sample]]
  xistHsM <- as.matrix(dgeDF[dgeDF$hgnc_symbol == "XIST", -c(1:3)])
  barplot(xistHsM
          , xaxt = "n", xlab = "Cells"
          , ylab = "Expression (counts)"
          , main = paste0(graphTitle
                          , "\n", sample, ": Drop-seq Xist expression"
                          , "\n", "> 10^4 human reads aligned, < 10^3 mouse reads aligned"))
}
dev.off()

## Log2 (expression + 1) by chromosome
for (sample in names(dgeLDF)) {
  dgeDF <- dgeLDF[[sample]]
  ggDF <- melt(dgeDF)
  # Remove biomart extra chromosomes
  rmvChr <- c(grep("CHR*", unique(ggDF$chromosome_name), value = TRUE)
              , grep("KI27*", unique(ggDF$chromosome_name), value = TRUE))
  ggDF <- ggDF[! ggDF$chromosome_name %in% rmvChr, ]
  ggDF$value <- log(ggDF$value + 1, 2)
  aggregate(ggDF$value, list(chr = ggDF$chromosome_name), mean)
  ggplot(ggDF, aes(y = value, x = chromosome_name)) +
    geom_jitter(alpha = 0.25) +
    ylab("Log2(counts + 1)") +
    xlab("Chromosome") + 
    ggtitle(paste0(graphTitle
                   , "\n", sample, ": Drop-seq expression by chromosome"
                   , "\n", "> 10^4 human reads aligned, < 10^3 mouse reads aligned"
                   , "\n"))
  ggsave(paste0(outGraph, "Expression_By_Chromosome_", sample, ".png"), width = 9, height = 6)
}
################################################################################

### Filter digital gene expression tables for > 20,000 human and < 1000 mouse reads

# Convert FASTQ line count to read count (divide by 4)
tReadsDF$TOTAL_FRAGMENTS <- tReadsDF$LINES_IN_FASTQ / 4

# Loop through sample directories, make list of filtered digital gene expression
# data frames
dgeLDF <- list()
for (dir in dirs) {
  
  # Load human counts pre UMI collapse from human mouse mapping to filter
  hsDat <- read.table(paste0(inParent, "/human_mouse/SxaQSEQsXap108L1/", dir
    , "/out_gene_exon_tagged_human.counts.summary.txt"), header = TRUE)
  
  print("")
  print(dir)
  str(hsDat)
  # List of barcodes with > X human exonic reads
  print(dim(hsDat[hsDat$NUM_TRANSCRIPTS > 20000, ]))
  hsBarcodes <- hsDat[hsDat$NUM_TRANSCRIPTS > 20000, ]$CELL_BARCODE
  
  # Load mouse counts pre UMI collapse from human mouse mapping to filter
  mmDat <- read.table(paste0(inParent, "/human_mouse/SxaQSEQsXap108L1/", dir
    , "/out_gene_exon_tagged_mouse.counts.summary.txt"), header = TRUE)
  str(mmDat)
  # List of barcodes with < X mouse exonic reads
  print(dim(mmDat[mmDat$NUM_TRANSCRIPTS < 1000, ]))
  mmBarcodes <- mmDat[mmDat$NUM_TRANSCRIPTS < 1000, ]$CELL_BARCODE
  
  hsDgeDF <- read.table(paste0(inParent
    , "/Homo_sapiens.GRCh37.75.dna.primary_assembly_NoERCC/SxaQSEQsXap108L1/"
    , dir, "/out_gene_exon_tagged_human.dge.txt.gz")
    , header = TRUE, row.names = 1)
  
  # Filtering digital gene expression table to list of cell barcodes with > X
  # human and < X mouse exonic reads
  hsFtDgeDF <- hsDgeDF[ ,colnames(hsDgeDF) %in% mmBarcodes & colnames(hsDgeDF) %in% hsBarcodes]
  print(dim(hsFtDgeDF))
  dgeLDF[[dir]] <- hsFtDgeDF
}
################################################################################

### Plot Xist and Y chromosome expression for >2e4 human reads

## Convert Ensembl IDs to Gene Symbols

for (sample in names(dgeLDF)) {
  dgeDF <- dgeLDF[[sample]]
  # Dataframe of ensembl IDs and chromosome
  # Human
  hsEnsemblChrDF <- AddChromosomeHuman(row.names(dgeDF))
  print(dim(hsEnsemblChrDF))
  
  # Add chromosome to expression dataframe
  # Human
  dgeLDF[[sample]] <- merge(hsEnsemblChrDF, dgeDF, by.x = "hgnc_symbol", by.y = "row.names")
}

## Y expression
pdf(paste0(outGraph, "Mean_Ychr_Expression_G2e4Hs.pdf"))
for (sample in names(dgeLDF)) {
  dgeDF <- dgeLDF[[sample]]
  yHsFtExDF <- dgeDF[dgeDF$chromosome_name == "Y", ]
  barplot(sort(apply(yHsFtExDF[ ,-c(1:3)], 2, mean))
    , xaxt = "n", xlab = "Cells"
    , ylab = "Expression (counts)"
    , main = paste0(graphTitle
      , "\n", sample, ": Drop-seq mean Y chromosome expression"
      , "\n", "> 2e4 human reads aligned, < 1e3 mouse reads aligned"))
}
dev.off()




