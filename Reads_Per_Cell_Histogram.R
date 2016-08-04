



dirs <- list.files("../data/bam/Homo_sapiens.GRCh37.75.dna.primary_assembly_NoERCC/SxaQSEQsXap108L1")

# dirs <- c("N702", "N703")

dir.create("../analysis/graphs", recursive = TRUE)

pdf(paste0("../analysis/graphs/Reads_Per_Cell_Histogram.pdf"))
for (dir in dirs) {
  a=read.table(paste0("../data/bam/Homo_sapiens.GRCh37.75.dna.primary_assembly_NoERCC/SxaQSEQsXap108L1/", dir
                      , "/out_cell_readcounts.txt.gz")
               , header=F, stringsAsFactors=F)
  
  x=cumsum(a$V1)
  x=x/max(x)
  plot(1:length(x), x, type='l', col="blue"
       , main = paste0("Reads Per Cell Histogram: ", dir)
       , xlab="Cell barcodes sorted by number of reads [descending]"
       , ylab="Cumulative fraction of reads", xlim=c(1, 5000))
}
dev.off()

