# Damon Polioudakis
# 2016-06-30
# Plot number of genes detected from Drop-seq digital gene expression
# Compare percent transcripts remaining after collapsing and not collapsing UMIs
################################################################################

rm(list=ls())
sessionInfo()

## Load data and assign variables

dirs <- c(paste0("human_mouse/SxaQSEQsXap108L1/"
                 , list.files("../data/digital_gene_expression/human_mouse/SxaQSEQsXap108L1"))
          , paste0("Homo_sapiens.GRCh37.75.dna.primary_assembly_NoERCC/SxaQSEQsXap108L1/"
                 , list.files("../data/digital_gene_expression/human_mouse/SxaQSEQsXap108L1")))
# dirs <- c("N702", "N703")

inParent <- "../data/digital_gene_expression"

dir.create("../analysis/graphs", recursive = TRUE)
################################################################################

## Digital Gene Expression

pdf("../analysis/graphs/Number_Genes_Transcripts_Detected_DGE_Genes.pdf")
for (dir in dirs) {
  dgeDF <- read.table(paste0(inParent, "/", dir
                          , "/out_gene_exon_tagged_human.dge.summary.txt")
                   , header = TRUE, stringsAsFactors=F)
  print(head(dgeDF))
  mn <- round(mean(dgeDF$NUM_GENES, na.rm = TRUE, ), 2)
  
  barplot(sort(dgeDF$NUM_GENES, decreasing = TRUE), col = "blue"
          # , names.arg = dgeDF$CELL_BARCODE
          , main = paste0("Number Genes Detected (Drop-seq Digital Expression): "
                          , "\n", dir
                          ,"\nMean: ", mn)
          , xlab = "Cell barcodes"
          , ylab = "Number of genes detected"
  )
}
dev.off()
################################################################################

## Counts Gene Expression

pdf("../analysis/graphs/Number_Genes_Transcripts_Detected_Counts_Genes.pdf")
for (dir in dirs) {          
  countDF <- read.table(paste0(inParent, "/", dir
                                     , "/out_gene_exon_tagged_human.counts.summary.txt")
                              , header = TRUE, stringsAsFactors=F)
  print(head(countDF))
  mn <- round(mean(countDF$NUM_GENES, na.rm = TRUE), 2)
          
  barplot(sort(countDF$NUM_GENES, decreasing = TRUE), col = "blue"
           # , names.arg = countDF$CELL_BARCODE
           , main = paste0("Number Genes Detected (Drop-seq Counts Expression): "
                                  , "\n", dir
                                  ,"\nMean: ", mn)
           , xlab = "Cell barcodes"
           , ylab = "Number of genes detected"
  )
}
dev.off()
################################################################################

## Digital Gene Expression Transcripts

pdf("../analysis/graphs/Number_Genes_Transcripts_Detected_DGE_Transcripts.pdf")
for (dir in dirs) {
  dgeDF <- read.table(paste0(inParent, "/", dir
                          , "/out_gene_exon_tagged_human.dge.summary.txt")
                   , header = TRUE, stringsAsFactors=F)
  print(head(dgeDF))
  mn <- round(mean(dgeDF$NUM_TRANSCRIPTS, na.rm = TRUE), 2)
  
  barplot(sort(dgeDF$NUM_TRANSCRIPTS, decreasing = TRUE), col = "blue"
          # , names.arg = dgeDF$CELL_BARCODE
          , main = paste0("Number Transcripts Detected (Drop-seq Digital Expression): "
                          , "\n", dir
                          ,"\nMean: ", mn)
          , xlab = "Cell barcodes"
          , ylab = "Number of transcripts detected"
  )
}
dev.off()
################################################################################

## Counts Expression Transcripts

pdf("../analysis/graphs/Number_Genes_Transcripts_Detected_Counts_Transcripts.pdf")
for (dir in dirs) {
  countDF <- read.table(paste0(inParent, "/", dir
                             , "/out_gene_exon_tagged_human.counts.summary.txt")
                      , header = TRUE, stringsAsFactors=F)
  print(head(countDF))
  mn <- round(mean(countDF$NUM_TRANSCRIPTS, na.rm = TRUE), 2)
  
  barplot(sort(countDF$NUM_TRANSCRIPTS, decreasing = TRUE), col = "blue"
          # , names.arg = countDF$CELL_BARCODE
          , main = paste0("Number Transcripts Detected (Drop-seq Counts Expression): "
                          , "\n", dir
                          ,"\nMean: ", mn)
          , xlab = "Cell barcodes"
          , ylab = "Number of transcripts detected"
  )
}
dev.off()
################################################################################

## Percent transcripts remaining after collapsing UMIs

pdf("../analysis/graphs/Number_Genes_Transcripts_Detected_Percent_Remaining_Transcripts.pdf")
for (dir in dirs) {
  dgeDF <- read.table(paste0(inParent, "/", dir
                             , "/out_gene_exon_tagged_human.dge.summary.txt")
                      , header = TRUE, stringsAsFactors=F)
  countDF <- read.table(paste0(inParent, "/", dir
                               , "/out_gene_exon_tagged_human.counts.summary.txt")
                        , header = TRUE, stringsAsFactors=F)
  
  df <- merge(dgeDF, countDF, by.x = "CELL_BARCODE", by.y = "CELL_BARCODE")
  df <- df[order(df$NUM_TRANSCRIPTS.x, decreasing = TRUE), ]
  pctRemain <- df$NUM_TRANSCRIPTS.x / df$NUM_TRANSCRIPTS.y * 100
  mn <- round(mean(pctRemain, na.rm = TRUE), 2)
  
  barplot(pctRemain, col = "blue"
          # , names.arg = countDF$CELL_BARCODE
          , main = paste0("Percent transcripts detected remaining after collapsing UMIs: "
                          , "\n", dir
                          ,"\nMean: ", mn)
          , xlab = "Cell barcodes"
          , ylab = "Percent remaining")
}
dev.off()
################################################################################
################################################################################

### Mouse cell barcodes subset

## Digital Gene Expression

pdf("../analysis/graphs/Number_Genes_Transcripts_Detected_DGE_Genes_MmSubset.pdf")
for (dir in dirs) {
  hsDat <- read.table(paste0(inParent, "/", dir, "/out_gene_exon_tagged_human.dge.summary.txt"), header = TRUE)
  
  str(hsDat)
  print(dim(hsDat[hsDat$NUM_TRANSCRIPTS < 250, ]))
  mmBarcodes <- hsDat[hsDat$NUM_TRANSCRIPTS < 250, ]$CELL_BARCODE
  
  dgeDF <- read.table(paste0(inParent, "/", dir
                          , "/out_gene_exon_tagged_mouse.dge.txt.gz")
                   , header = TRUE, stringsAsFactors=F)
                   
  idxMm <- colnames(dgeDF) %in% mmBarcodes
  dgeMmDF <- dgeDF[ ,idxMm]
  dgeMmDF <- dgeMmDF >= 1
  # Convert to numeric
  dgeMmDF <- dgeMmDF + 0
  nGenes <- colSums(dgeMmDF)

  mn <- round(mean(nGenes, na.rm = TRUE, ), 2)
  
  barplot(sort(nGenes, decreasing = TRUE), col = "blue"
          # , names.arg = dgeDF$CELL_BARCODE
          , main = paste0("Number Genes Detected (Drop-seq Digital Expression): "
          				  , "\nSubset to mouse cells (<250 human transcripts)"
                          , "\n", dir
                          , "\nMean: ", mn)
          , xlab = "Cell barcodes"
          , ylab = "Number of genes detected"
  )
}
dev.off()
################################################################################

## Counts Gene Expression

pdf("../analysis/graphs/Number_Genes_Transcripts_Detected_Counts_Genes_MmSubset.pdf")
for (dir in dirs) {
  hsDat <- read.table(paste0(inParent, "/", dir, "/out_gene_exon_tagged_human.dge.summary.txt"), header = TRUE)
  
  str(hsDat)
  print(dim(hsDat[hsDat$NUM_TRANSCRIPTS < 250, ]))
  mmBarcodes <- hsDat[hsDat$NUM_TRANSCRIPTS < 250, ]$CELL_BARCODE
  
  dgeDF <- read.table(paste0(inParent, "/", dir
                          , "/out_gene_exon_tagged_mouse.counts.txt.gz")
                   , header = TRUE, stringsAsFactors=F)
                   
  idxMm <- colnames(dgeDF) %in% mmBarcodes
  dgeMmDF <- dgeDF[ ,idxMm]
  dgeMmDF <- dgeMmDF >= 1
  # Convert to numeric
  dgeMmDF <- dgeMmDF + 0
  nGenes <- colSums(dgeMmDF)

  mn <- round(mean(nGenes, na.rm = TRUE, ), 2)
  
  barplot(sort(nGenes, decreasing = TRUE), col = "blue"
          # , names.arg = dgeDF$CELL_BARCODE
          , main = paste0("Number Genes Detected (Drop-seq Counts Expression): "
          				  , "\nSubset to mouse cells (<250 human transcripts)"
                          , "\n", dir
                          , "\nMean: ", mn)
          , xlab = "Cell barcodes"
          , ylab = "Number of genes detected"
  )
}
dev.off()

# # for (dir in dirs) {          
  # countDF <- read.table(paste0("../data/digital_gene_expression/human_mouse/SxaQSEQsXap108L1/", dir
                                     # , "/out_gene_exon_tagged_human.counts.summary.txt")
                              # , header = TRUE, stringsAsFactors=F)
  # print(head(countDF))
  # mn <- round(mean(countDF$NUM_GENES, na.rm = TRUE), 2)
          
  # barplot(sort(countDF$NUM_GENES, decreasing = TRUE), col = "blue"
           # # , names.arg = countDF$CELL_BARCODE
           # , main = paste0("Number Genes Detected (Drop-seq Counts Expression): "
                                  # , dir
                                  # ,"\nMean: ", mn)
           # , xlab = "Cell barcodes"
           # , ylab = "Number of genes detected"
  # )
# }
# dev.off()
################################################################################

## Digital Gene Expression Transcripts

pdf("../analysis/graphs/Number_Genes_Transcripts_Detected_DGE_Transcripts_MmSubset.pdf")
for (dir in dirs) {
  hsDat <- read.table(paste0(inParent, "/", dir, "/out_gene_exon_tagged_human.dge.summary.txt"), header = TRUE)
  
  str(hsDat)
  print(dim(hsDat[hsDat$NUM_TRANSCRIPTS < 250, ]))
  mmBarcodes <- hsDat[hsDat$NUM_TRANSCRIPTS < 250, ]$CELL_BARCODE
  
  dgeDF <- read.table(paste0(inParent, "/", dir
                          , "/out_gene_exon_tagged_mouse.dge.txt.gz")
                   , header = TRUE, stringsAsFactors=F)
                   
  idxMm <- colnames(dgeDF) %in% mmBarcodes
  dgeMmDF <- dgeDF[ ,idxMm]
  nTrscpts <- colSums(dgeMmDF)

  mn <- round(mean(nTrscpts, na.rm = TRUE, ), 2)
  
  barplot(sort(nTrscpts, decreasing = TRUE), col = "blue"
          # , names.arg = dgeDF$CELL_BARCODE
          , main = paste0("Number Transcripts Detected (Drop-seq Digital Expression): "
          				  , "\nSubset to mouse cells (<250 human transcripts)"
                          , "\n", dir
                          , "\nMean: ", mn)
          , xlab = "Cell barcodes"
          , ylab = "Number of transcripts detected"
  )
}
dev.off()

# for (dir in dirs) {
  # dgeDF <- read.table(paste0("../data/digital_gene_expression/human_mouse/SxaQSEQsXap108L1/", dir
                          # , "/out_gene_exon_tagged_human.dge.summary.txt")
                   # , header = TRUE, stringsAsFactors=F)
  # print(head(dgeDF))
  # mn <- round(mean(dgeDF$NUM_TRANSCRIPTS, na.rm = TRUE), 2)
  
  # barplot(sort(dgeDF$NUM_TRANSCRIPTS, decreasing = TRUE), col = "blue"
          # # , names.arg = dgeDF$CELL_BARCODE
          # , main = paste0("Number Transcripts Detected (Drop-seq Digital Expression): "
                          # , dir
                          # ,"\nMean: ", mn)
          # , xlab = "Cell barcodes"
          # , ylab = "Number of transcripts detected"
  # )
# }
# dev.off()
################################################################################

## Counts Expression Transcripts

pdf("../analysis/graphs/Number_Genes_Transcripts_Detected_Counts_Transcripts_MmSubset.pdf")
for (dir in dirs) {
  hsDat <- read.table(paste0(inParent, "/", dir, "/out_gene_exon_tagged_human.dge.summary.txt"), header = TRUE)
  
  str(hsDat)
  print(dim(hsDat[hsDat$NUM_TRANSCRIPTS < 250, ]))
  mmBarcodes <- hsDat[hsDat$NUM_TRANSCRIPTS < 250, ]$CELL_BARCODE
  
  dgeDF <- read.table(paste0(inParent, "/", dir
                          , "/out_gene_exon_tagged_mouse.counts.txt.gz")
                   , header = TRUE, stringsAsFactors=F)
                   
  idxMm <- colnames(dgeDF) %in% mmBarcodes
  dgeMmDF <- dgeDF[ ,idxMm]
  nTrscpts <- colSums(dgeMmDF)

  mn <- round(mean(nTrscpts, na.rm = TRUE, ), 2)
  
  barplot(sort(nTrscpts, decreasing = TRUE), col = "blue"
          # , names.arg = dgeDF$CELL_BARCODE
          , main = paste0("Number Transcripts Detected (Drop-seq Counts Expression): "
          				  , "\nSubset to mouse cells (<250 human transcripts)"
                          , "\n", dir
                          , "\nMean: ", mn)
          , xlab = "Cell barcodes"
          , ylab = "Number of transcripts detected"
  )
}
dev.off()

# for (dir in dirs) {
  # countDF <- read.table(paste0("../data/digital_gene_expression/human_mouse/SxaQSEQsXap108L1/", dir
                             # , "/out_gene_exon_tagged_human.counts.summary.txt")
                      # , header = TRUE, stringsAsFactors=F)
  # print(head(countDF))
  # mn <- round(mean(countDF$NUM_TRANSCRIPTS, na.rm = TRUE), 2)
  
  # barplot(sort(countDF$NUM_TRANSCRIPTS, decreasing = TRUE), col = "blue"
          # # , names.arg = countDF$CELL_BARCODE
          # , main = paste0("Number Transcripts Detected (Drop-seq Counts Expression): "
                          # , dir
                          # ,"\nMean: ", mn)
          # , xlab = "Cell barcodes"
          # , ylab = "Number of transcripts detected"
  # )
# }
# dev.off()
################################################################################

## Percent transcripts remaining after collapsing UMIs

pdf("../analysis/graphs/Number_Genes_Transcripts_Detected_Percent_Remaining_Transcripts_MmSubset.pdf")
for (dir in dirs) {
	
  hsDat <- read.table(paste0(inParent, "/", dir, "/out_gene_exon_tagged_human.dge.summary.txt"), header = TRUE)
  
  str(hsDat)
  print(dim(hsDat[hsDat$NUM_TRANSCRIPTS < 250, ]))
  mmBarcodes <- hsDat[hsDat$NUM_TRANSCRIPTS < 250, ]$CELL_BARCODE

  dgeDF <- read.table(paste0(inParent, "/", dir
                             , "/out_gene_exon_tagged_mouse.dge.txt.gz")
                      , header = TRUE, stringsAsFactors=F)
  countDF <- read.table(paste0(inParent, "/", dir
                               , "/out_gene_exon_tagged_mouse.counts.txt.gz")
                        , header = TRUE, stringsAsFactors=F)
                        
  idxMm <- colnames(dgeDF) %in% mmBarcodes
  dgeMmDF <- dgeDF[ ,idxMm]
  nTrscptsDGE <- data.frame(colSums(dgeMmDF))
  
  idxMm <- colnames(countDF) %in% mmBarcodes
  countMmDF <- countDF[ ,idxMm]
  nTrscptsCount <- data.frame(colSums(countMmDF))
  
  print(nTrscptsCount)
  
  df <- merge(nTrscptsDGE, nTrscptsCount, by.x = "row.names", by.y = "row.names")
  
  print(str(df))

  df <- df[order(df$colSums.dgeMmDF., decreasing = TRUE), ]
  pctRemain <- df$colSums.dgeMmDF. / df$colSums.countMmDF. * 100
  mn <- round(mean(pctRemain, na.rm = TRUE), 2)
  
  barplot(pctRemain, col = "blue"
          # , names.arg = countDF$CELL_BARCODE
          , main = paste0("Percent transcripts detected remaining after collapsing UMIs: "
                          , "\n", dir
                          ,"\nMean: ", mn)
          , xlab = "Cell barcodes"
          , ylab = "Percent remaining")
}
dev.off()
################################################################################
################################################################################

### Subset to human cells (Mouse < 1000 & Human > 10,000)

## Digital Gene Expression

pdf("../analysis/graphs/Number_Genes_Transcripts_Detected_DGE_Genes_SubsetHs10e4Mm1000.pdf")
for (dir in sub(pattern = "^(.*?)/", replacement = "", x = dirs[1:3])) {
  hsDat <- read.table(paste0(inParent, "/human_mouse/", dir
                             , "/out_gene_exon_tagged_human.counts.summary.txt"), header = TRUE)
  
  print("")
  print(dir)
  str(hsDat)
  print(dim(hsDat[hsDat$NUM_TRANSCRIPTS > 10000, ]))
  hsBarcodes <- hsDat[hsDat$NUM_TRANSCRIPTS > 10000, ]$CELL_BARCODE
  
  mmDat <- read.table(paste0(inParent, "/human_mouse/", dir
                             , "/out_gene_exon_tagged_mouse.counts.summary.txt"), header = TRUE)
  str(mmDat)
  print(dim(mmDat[mmDat$NUM_TRANSCRIPTS < 1000, ]))
  mmBarcodes <- mmDat[mmDat$NUM_TRANSCRIPTS < 1000, ]$CELL_BARCODE
  
  dgeDF <- read.table(paste0(inParent, "/Homo_sapiens.GRCh37.75.dna.primary_assembly_NoERCC/", dir
                             , "/out_gene_exon_tagged_human.dge.txt.gz")
                      , header = TRUE, stringsAsFactors=F)
  
  idxHs <- ! colnames(dgeDF) %in% mmBarcodes
  idxHs <- colnames(dgeDF) %in% hsBarcodes
  print(table(idxHs))
  dgeHsDF <- dgeDF[ ,idxHs]
  dgeHsDF <- dgeHsDF >= 1
  # Convert to numeric
  dgeHsDF <- dgeHsDF + 0
  nGenes <- colSums(dgeHsDF)
  
  mn <- round(mean(nGenes, na.rm = TRUE), 2)
  
  barplot(sort(nGenes, decreasing = TRUE), col = "blue"
          # , names.arg = dgeDF$CELL_BARCODE
          , main = paste0("Number Genes Detected (Drop-seq Digital Expression): "
                          , "\nSubset to human cells (>10,000 human trancripts and <1000 mouse transcripts)"
                          , "\n", dir
                          , "\nMean: ", mn)
          , xlab = "Cell barcodes"
          , ylab = "Number of genes detected"
  )
}
dev.off()
################################################################################



## Code for testing

# dgeStatsDF <- read.table("../data/digital_gene_expression/human_mouse/SxaQSEQsXap108L1/N701/out_gene_exon_tagged_human.dge.summary.txt"
#                          , header = TRUE)
# 
# # dgeDatDF <- read.table("../data/digital_gene_expression/human_mouse/SxaQSEQsXap108L1/N701/out_gene_exon_tagged_human.dge.txt.gz"
# # , header = TRUE)
# 
# countsStatsDF <- read.table("../data/digital_gene_expression/human_mouse/SxaQSEQsXap108L1/N701/out_gene_exon_tagged_human.counts.summary.txt"
#                             , header = TRUE)
# 
# # countsDatDF <- read.table("../data/digital_gene_expression/human_mouse/SxaQSEQsXap108L1/N701/out_gene_exon_tagged_human.counts.txt.gz"
# # , header = TRUE)
# 
# df <- merge(dgeStatsDF, countsStatsDF, by.x = "CELL_BARCODE", by.y = "CELL_BARCODE")
# df <- df[order(df$NUM_TRANSCRIPTS.x, decreasing = TRUE), ]
# pctRemain <- df$NUM_TRANSCRIPTS.x / df$NUM_TRANSCRIPTS.y * 100

# sort(colSums(dgeDatDF[ ,-1]))
