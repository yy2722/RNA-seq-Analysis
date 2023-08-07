#!/bin/bash
#$ -j n
#$ -N olego_human
#$ -t 1-23
#$ -tc 5
#$ -o /mnt/chromatin/home/yy2722/logs/olegoMapping
#$ -e /mnt/chromatin/home/yy2722/logs/olegoMapping

source /mnt/chromatin/home/yy2722/.bash_profile
pwd
conda activate quantas

sraAcc=$(cat Yoshimi_humanData.txt | awk "NR==$SGE_TASK_ID")

homePath="/mnt/chromatin/home/yy2722/"

# directories for reference genome files
genome_path="${homePath}data/genome/hg19/"
genome_cfg="${genome_path}hg.cfg"
genome_hmrBed="${genome_path}hg19.intron.hmr.bed"
genome_olego="${genome_path}hg19_olego"

# olego input and output directories
olego_oF="${homePath}olegoOutput/${sraAcc}f.sam"
olego_iF="${homePath}fastqFiles/${sraAcc}_1.fastq"
olego_oR="${homePath}olegoOutput/${sraAcc}r.sam"
olego_iR="${homePath}fastqFiles/${sraAcc}_2.fastq"

# map forward strand reads to the reference genome
olego -v -t 16 -r $genome_cfg -j $genome_hmrBed -o $olego_oF $genome_olego $olego_iF

# map reverse strand reads
olego -v -t 16 -r $genome_cfg -j $genome_hmrBed -o $olego_oR $genome_olego $olego_iR
