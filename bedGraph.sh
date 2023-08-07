#!/bin/bash
#$ -j n
#$ -N bedGraph_mouse
#$ -t 1-6
#$ -tc 10
#$ -o /mnt/chromatin/home/yy2722/logs/bedGraph
#$ -e /mnt/chromatin/home/yy2722/logs/bedGraph

source /mnt/chromatin/home/yy2722/.bash_profile
pwd
conda activate quantas

sraAcc=$(cat Yoshimi_mouseData.txt | awk "NR==$SGE_TASK_ID")

homePath="/mnt/chromatin/home/yy2722/quantasOutput/"

tag2profile.pl -v -big -weight -c ${homePath}${sraAcc}_cache -exact -of bedgraph -n "${sraAcc}_bed_graph" ${homePath}${sraAcc}_gapless_out/pair.gapless.bed ${homePath}${sraAcc}.bedGraph
