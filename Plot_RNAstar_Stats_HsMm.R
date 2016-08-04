



starDF <- read.table("../metadata/human_mouse/RNAstar_Stats_SxaQSEQsXap108L1.txt"
           , header = TRUE, fill = TRUE, sep = "\t")

str(starDF)

df <- data.frame(t(starDF[ ,c(1, 2, 3, 5, 6, 7, 20, 21, 22, 23, 25, 26, 27)]))

rNames <- gsub("X", "%", row.names(df))
rNames <- gsub("\\.", " ", rNames)

outDF <- cbind(rNames, df)

dir.create("../analysis/tables", recursive = TRUE)
write.table(outDF, "../analysis/tables/Plot_RNAstar_Stats_HsMm.txt", quote = FALSE
            , sep = "\t", row.names = FALSE, col.names = FALSE)
