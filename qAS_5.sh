#!/bin/bash
#$ -j n
#$ -N quantifyAS_human
#$ -t 1-23
#$ -tc 5
#$ -o /mnt/chromatin/home/yy2722/logs/quantifyAS
#$ -e /mnt/chromatin/home/yy2722/logs/quantifyAS

source /mnt/chromatin/home/yy2722/.bash_profile
pwd
conda activate quantas

sraAcc=$(cat Yoshimi_humanData.txt | awk "NR==$SGE_TASK_ID")

homePath="/mnt/chromatin/home/yy2722/"
cachePath="${homePath}quantasOutput/${sraAcc}_cache"

# directories for reference genome files
genome_path="${homePath}data/genome/hg19/"

genome_conf="${genome_path}hg19.conf"
geneSymbol="${genome_path}hg19.exon.uniq.core.id2gene2symbol"

ASinput="${homePath}quantasOutput/${sraAcc}_gapless_out/pair.gapless.bed"
ASoutput="${homePath}quantasOutput/${sraAcc}_countit_out"

summarize_splicing_wrapper.pl -c $cachePath -v -big -weight -conf $genome_conf -dbkey hg19 -cass -taca -alt5 -alt3 -mutx -iret $ASinput $ASoutput
