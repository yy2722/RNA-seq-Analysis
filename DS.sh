#!/bin/bash
#$ -j n
#$ -N splicingDiff_human
#$ -t 1-6
#$ -e /mnt/chromatin/home/yy2722/logs/splicingDiff
#$ -o /mnt/chromatin/home/yy2722/logs/splicingDiff

source /mnt/chromatin/home/yy2722/.bash_profile
pwd
conda activate quantas

AStype=("cass" "taca" "alt5" "alt3" "mutx" "iret")

homeDir="/mnt/chromatin/home/yy2722/"

test_splicing_diff.pl -type ${AStype[$SGE_TASK_ID-1]} -v --min-cov 20 --id2gene2symbol ${homeDir}data/genome/hg19/annotation/Hs.seq.all.AS.chrom.can.id2gene2symbol ${homeDir}quantasOutput/hs.dataset.group.conf ${homeDir}quantasOutput/hs.${AStype[$SGE_TASK_ID-1]}_dataset.diff.txt
