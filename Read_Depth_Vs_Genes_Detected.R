




dirs <- list.files("../data/bam/human_mouse/SxaQSEQsXap108L1")
# dir <- c("N702")

dir.create("../analysis/graphs", recursive = TRUE)

pdf("../analysis/graphs/Read_Depth_Vs_Genes_Detected.pdf")
for (dir in dirs) {
  rdDF <- read.table(paste0("../data/bam/human_mouse/SxaQSEQsXap108L1/", dir
                          , "/out_cell_readcounts.txt.gz")
                   , header = FALSE, stringsAsFactors = F, nrows = 400)
  colnames(rdDF) <- c("READ_DEPTH", "CELL_BARCODE")
  print(head(rdDF))
  
  dgeDF <- read.table(paste0("../data/digital_gene_expression/human_mouse/SxaQSEQsXap108L1/", dir
                             , "/out_gene_exon_tagged_human.dge.summary.txt")
                      , header = TRUE, stringsAsFactors=F)
  print(head(dgeDF))
  
  df <- merge(rdDF, dgeDF, by.x = "CELL_BARCODE", by.y = "CELL_BARCODE")
  print(head(df))
  df <- df[order(df$READ_DEPTH), ]
  
  plot(df$READ_DEPTH, df$NUM_GENES
       , main = paste0("Read_Depth_Vs_Genes_Detected.R"
                         , "\n", dir)
       , xlab = "Read Depth"
       , ylab = "Number of genes detected"
       )
  plot(df$READ_DEPTH, df$NUM_GENES
       , main = paste0("Read_Depth_Vs_Genes_Detected.R"
                       , "\nX-limit adjusted"
                       , "\n", dir)
       , xlim = c(0, 10^5)
       , xlab = "Read Depth"
       , ylab = "Number of genes detected"
  )
  
}
dev.off()