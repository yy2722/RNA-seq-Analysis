#!/bin/bash
#$ -j n
#$ -N fasterq_human
#$ -t 4-23
#$ -tc 5
#$ -o /mnt/chromatin/home/yy2722/logs/fasterq
#$ -e /mnt/chromatin/home/yy2722/logs/fasterq

source /mnt/chromatin/home/yy2722/.bash_profile
pwd
conda activate quantas

# get SRA number from txt file
sraAcc=$(cat Yoshimi_humanData.txt | awk "NR==$SGE_TASK_ID")

fasterq-dump -O /mnt/chromatin/home/yy2722/fastqFiles  $sraAcc
