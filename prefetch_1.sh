#!/bin/bash
#$ -j n
#$ -N prefetch_human
#$ -t 4-23
#$ -tc 8
#$ -o /mnt/chromatin/home/yy2722/logs/prefetch
#$ -e /mnt/chromatin/home/yy2722/logs/prefetch

source /mnt/chromatin/home/yy2722/.bash_profile
pwd
conda activate quantas

# get SRA number from txt file
sraAcc=$(cat Yoshimi_humanData.txt | awk "NR==$SGE_TASK_ID")

# run prefetch
prefetch $sraAcc

# validate data integrity
vdb-validate /mnt/chromatin/home/yy2722/ncbi_download/sra/${sraAcc}.sra
