# Damon Polioudakis
# 2016-06-03
# Process RNA-seq data from QSEQ to demultplexed fastq and remove reads that are not PF=1
# Adopted form Jason Stein's Preprocess.R script
# Use QSUB .sh script to run
# Requires:
# 	0.1_Pair.bash 
# 	0.2_fastq.R1.qual.convert.pl 
# 	0.3_fastq.R2.qual.convert.pl 

###################################

##Load libraries
library(xlsx)
options(stringsAsFactors = FALSE)

###################################
##User set variables
#args = commandArgs(trailingOnly = TRUE)
##Sample information spreadsheet
sampleinfofname = "../metadata/Nextera_Indexes.xlsx"
#sampleinfofname = args[1]
cat('Sample info file:', sampleinfofname, '\n')
##Directory which contains the pools
pooldir = "../data/qseq"
##pooldir = args[2]
cat('Pool directory:', pooldir, '\n')

###################################
##Directory structure
##Directory which contains the qseq files
qseqdir = list.files(pooldir, pattern = "^SxaQSEQs[[:alnum:]]{6}L[[:digit:]]{1}$")
lane = substr(qseqdir, nchar(qseqdir[1]) - 1, nchar(qseqdir[1]))
qseqdir = list.files(pooldir, pattern = "^SxaQSEQs[[:alnum:]]{6}L[[:digit:]]{1}$", full.names = TRUE)
qseqLaneDir = list.files(pooldir, pattern = "^SxaQSEQs[[:alnum:]]{6}L[[:digit:]]{1}$")

###################################
##Program variables
nlanes = length(qseqdir)
cat('Number of lanes:', nlanes, '\n')

##Loop over lanes (if multiple lanes exist)
for (l in 1:nlanes) {
	
	cat("\nLane:", l)

	###################################
	##Read in spreadsheet containing sample information
	sampleinfo = read.xlsx(sampleinfofname, 1)
	###################################
	
    cat("\nWorking on directory:", qseqdir[l], "\n")
    outdir = paste0('../data/fastq/', qseqLaneDir[l])
    cat("Out directory:", outdir, "\n")
    dir.create(outdir, showWarnings = FALSE)

    ###################################
    # ##Merge two paired reads from qseq files 
    # cat('Merging two paired reads from qseq files\n')
    # dir.create(outdir, recursive = TRUE)
    # system(paste0("/geschwindlabshares/RNAseq_singlecellfetal/DS-001_DP/code/0.1_Pair.bash ", qseqdir[l], " ", outdir))
    
    # ##merge all qseq files into 1 file
    # cat('Merging all qseq files into one large file\n')
    # system(paste0("cat ", outdir, "/s* > ", outdir, "/merged.txt"))
    
    # ##Determine the number of reads in the lane
    # cat('Counting the number of reads in the lane\n')
    # unfilteredreads = system(paste0("wc -l ", outdir, "/merged.txt"), intern = TRUE)
    # unfilteredreads = as.numeric(unlist(strsplit(unfilteredreads, " "))[1])
    # cat('Number of unfiltered reads in the lane:', unfilteredreads, '\n')
    
    # ##Filter reads with PF = 1 (QC step, see manual from core)
    # cat('Filtering the reads to take only PF = 1\n')
    # # awk '$7 designates PF column in intermediate file after merging qseq files
    # system(paste0("cat ", outdir, "/merged.txt | awk '$7 == 1' > ", outdir, "/merged.PF1.txt"))
    # filteredreads = system(paste0("wc -l ", outdir, "/merged.PF1.txt"), intern = TRUE)
    # filteredreads = as.numeric(unlist(strsplit(filteredreads, " "))[1])
    # cat('Number of filtered reads in the lane:', filteredreads, '\n')
    # cat('Percentage of useable reads:', signif(filteredreads/unfilteredreads, 2)*100, '%\n')

    ###################################
    ##Loop over all samples within a pool
    for (i in 1:length(sampleinfo$Index_Name)) {

        ###################################
        ##De-multiplex data by indexes
        sampleIndexI7 <- gsub("TruSeq_[[:digit:]]*_", "", sampleinfo$i7_Index_Sequence[i])
        cat("QSUBing Sample", sampleinfo$Index_Name[i], "for de-muliplexing and converting to FASTQ...\n")
        cat(paste0("grep -P \'^", sampleIndexI7, "\\t\' ", outdir, "/merged.PF1.txt | cut -f  1-10 > ", outdir, "/", sampleinfo$Index_Name[i], ".txt\n"), file = paste0(outdir, "/Make_fastq_Sample_", sampleinfo$Index_Name[i], ".sh"))
		        
       ###################################
		##Convert to FASTQ and gzip the result
		cat(paste0("/geschwindlabshares/RNAseq_singlecellfetal/DS-001_DP/code/0.2_fastq.R1.qual.convert.pl ", outdir, "/", sampleinfo$Index_Name[i], ".txt > ", outdir, "/", sampleinfo$Index_Name[i], ".R1.fastq\n"), file = paste0(outdir, "/Make_fastq_Sample_", sampleinfo$Index_Name[i], ".sh"), append = TRUE)
		cat(paste0("/geschwindlabshares/RNAseq_singlecellfetal/DS-001_DP/code/0.3_fastq.R2.qual.convert.pl ", outdir, "/", sampleinfo$Index_Name[i], ".txt > ", outdir, "/", sampleinfo$Index_Name[i], ".R2.fastq\n"), file = paste0(outdir, "/Make_fastq_Sample_", sampleinfo$Index_Name[i], ".sh"), append = TRUE)
		cat(paste0("gzip -f ", outdir, "/", sampleinfo$Index_Name[i], ".R1.fastq\n"), file = paste0(outdir, "/Make_fastq_Sample_", sampleinfo$Index_Name[i], ".sh"), append = TRUE)
		cat(paste0("gzip -f ", outdir, "/", sampleinfo$Index_Name[i], ".R2.fastq\n"), file = paste0(outdir, "/Make_fastq_Sample_", sampleinfo$Index_Name[i], ".sh"), append = TRUE)
	
	# Deleting intermediate files with script may mess up fastq qsub - delete manually after
	# ##Remove all intermediate files
	# if (i == length(sampleinfo$Index_Name)) {
           # ###################################
           # ##Remove all intermediate files
           # # N*.txt is from Index_Name in indexes excel file (e.g. N718)
           # cat(paste0("rm ", outdir, "/N*.txt\n"), file = paste0(outdir, "/Make_fastq_Sample_", sampleinfo$Index_Name[i], ".sh"), append = TRUE)
           # cat(paste0("rm ", outdir, "/merged.txt\n"), file = paste0(outdir, "/Make_fastq_Sample_", sampleinfo$Index_Name[i], ".sh"), append = TRUE)
           # cat(paste0("rm ", outdir, "/merged.PF1.txt\n"), file = paste0(outdir, "/Make_fastq_Sample_", sampleinfo$Index_Name[i], ".sh"), append = TRUE)
        # }
	
	##Run this job through qsub
	system(paste0("chmod 775 ", outdir, "/Make_fastq_Sample_", sampleinfo$Index_Name[i], ".sh"))
	system(paste0("qsub -S /bin/bash -N SC_", sampleinfo$Index_Name[i], " -o ", outdir, "/Make_fastq_Sample_", sampleinfo$Index_Name[i], ".log -e ", outdir, "/Make_fastq_Sample_", sampleinfo$Index_Name[i], ".error -pe shared 8 -j y -cwd -m n -q geschwind.q ", outdir, "/Make_fastq_Sample_", sampleinfo$Index_Name[i], ".sh"), intern = TRUE)
	#jobids[i, l] = unlist(strsplit(jobids[i, l], " ", ))[3]
    }
}