# Script to submit an Rscript to Orion queue
# Use this qsub:
# qsub -cwd -V -N QtoF -S /bin/bash -l h_data=16G,h_rt=24:00:00 -o logs/0_QSEQ_to_fastq_QSUB_$(date +%Y%m%d)_demultiplex.log -e logs/0_QSEQ_to_fastq_QSUB_$(date +%Y%m%d)_demultiplex.error 0_QSEQ_to_fastq_QSUB.sh

#!/bin/bash

# Path to R-3.2.2 to use xlsx package
/share/apps/R-3.2.2/bin/Rscript --vanilla 0_QSEQ_to_fastq.R
