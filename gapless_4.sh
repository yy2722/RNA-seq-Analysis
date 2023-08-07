#!/bin/bash
#$ -j n
#$ -N gapless_human
#$ -t 1-23
#$ -tc 5
#$ -o /mnt/chromatin/home/yy2722/logs/gapless
#$ -e /mnt/chromatin/home/yy2722/logs/gapless

source /mnt/chromatin/home/yy2722/.bash_profile
pwd
conda activate quantas

sraAcc=$(cat Yoshimi_humanData.txt | awk "NR==$SGE_TASK_ID")

homePath="/mnt/chromatin/home/yy2722/"

# directories for reference genome files
genome_path="${homePath}data/genome/hg19/"
genome_exonTrio="${genome_path}hg19.exon.trio.hmr.nr.bed"

# olego directories
olegoF="${homePath}olegoOutput/${sraAcc}f.sam"
olegoR="${homePath}olegoOutput/${sraAcc}r.sam"

gapless_out="${homePath}quantasOutput/${sraAcc}_gapless_out"

# infer transcript structure
gapless_huge_file.pl -v -sam -uniq --split-size 10000000 -isoform $genome_exonTrio -E 400 -big --print-singleton -o $gapless_out $olegoF $olegoR  
