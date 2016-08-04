# Damon Polioudakis
# 2016-06-30
# Plot barnyard plots for both transcripts detected and genes detected
# Transcripts detected is reads mapping to exons after UMI collapse
################################################################################

rm(list=ls())
sessionInfo()

inDir <- "../data/digital_gene_expression/human_mouse/SxaQSEQsXap108L1"

dir.create("../analysis/graphs", recursive = TRUE)
################################################################################

## Genes detected

dirs <- list.files("../data/bam/human_mouse/SxaQSEQsXap108L1")
# dirs <- c("N702")

pdf("../analysis/graphs/Barnyard_Plots_Genes.pdf")
for (dir in dirs) {
	hsDat <- read.table(paste0(inDir, "/", dir, "/out_gene_exon_tagged_human.dge.summary.txt"), header = TRUE)
	mmDat <- read.table(paste0(inDir, "/", dir, "/out_gene_exon_tagged_mouse.dge.summary.txt"), header = TRUE)
	hsMmDat <- merge(hsDat, mmDat, by.x = "CELL_BARCODE", by.y = "CELL_BARCODE")
	colnames(hsMmDat) <- c("CELL_BARCODE", "NUM_GENES_Hs", "NUM_TRANSCRIPTS_Hs", "NUM_GENES_Mm", "NUM_TRANSCRIPTS_Mm")
	
	plot(hsMmDat$NUM_GENES_Hs, hsMmDat$NUM_GENES_Mm
		, main = paste0("Barnyard_Plot.R"
                         , "\nNumber of genes detected (digital gene expression)"
                         , "\n", dir)
       , xlab = "Human number of genes detected"
       , ylab = "Mouse number of genes detected"
       )
}
dev.off()


## Transcripts detected

pdf("../analysis/graphs/Barnyard_Plots_Transcripts.pdf")
for (dir in dirs) {
	hsDat <- read.table(paste0(inDir, "/", dir, "/out_gene_exon_tagged_human.dge.summary.txt"), header = TRUE)
	mmDat <- read.table(paste0(inDir, "/", dir, "/out_gene_exon_tagged_mouse.dge.summary.txt"), header = TRUE)
	hsMmDat <- merge(hsDat, mmDat, by.x = "CELL_BARCODE", by.y = "CELL_BARCODE")
	colnames(hsMmDat) <- c("CELL_BARCODE", "NUM_GENES_Hs", "NUM_TRANSCRIPTS_Hs", "NUM_GENES_Mm", "NUM_TRANSCRIPTS_Mm")

	plot(hsMmDat$NUM_TRANSCRIPTS_Hs, hsMmDat$NUM_TRANSCRIPTS_Mm
		, main = paste0("Barnyard_Plot.R"
                         , "\nNumber of transcripts detected (digital gene expression)"
                         , "\n", dir)
       , xlab = "Human number of transcripts detected"
       , ylab = "Mouse number of transcripts detected"
       )
    # Line added for subsetting mouse cells
    plot(hsMmDat$NUM_TRANSCRIPTS_Hs, hsMmDat$NUM_TRANSCRIPTS_Mm
    	, abline(v = 250, col = 2)
		, main = paste0("Barnyard_Plot.R"
                         , "\nNumber of transcripts detected (digital gene expression)"
                         , "\n", dir)
       , xlab = "Human number of transcripts detected"
       , ylab = "Mouse number of transcripts detected"
       )
    )
}
dev.off()