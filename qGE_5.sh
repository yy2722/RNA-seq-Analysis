#!/bin/bash
#$ -j n
#$ -N DE_human
#$ -t 1-23
#$ -tc 5
#$ -o /mnt/chromatin/home/yy2722/logs/quantifyDE
#$ -e /mnt/chromatin/home/yy2722/logs/quantifyDE

source /mnt/chromatin/home/yy2722/.bash_profile
pwd
conda activate quantas

sraAcc=$(cat Yoshimi_humanData.txt | awk "NR==$SGE_TASK_ID")

homePath="/mnt/chromatin/home/yy2722/"

# directories for reference genome files
genome_path="${homePath}data/genome/hg19/"

coreBed="${genome_path}hg19.exon.uniq.core.bed"
geneSymbol="${genome_path}hg19.exon.uniq.core.id2gene2symbol"

# output directory
DEinput="${homePath}quantasOutput/${sraAcc}_gapless_out/pair.gapless.bed"
DEoutput="${homePath}quantasOutput/${sraAcc}.expr.txt"

summarize_expression_wrapper.pl -big --cache ${homePath}quantasOutput/${sraAcc}_cache -exon $coreBed -e2g $geneSymbol -weight -v $DEinput $DEoutput
