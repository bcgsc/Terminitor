#!/bin/bash
#SBATCH --ntasks=48
#SBATCH -N 1-1
#SBATCH --exclude='n[121-124,129-132,224-228]'
#SBATCH --mem=360G

sample=$1
/usr/bin/time -v minimap2 -ax splice hg38.mmi -t 12 ../assembly/${sample}_combined/rnabloom.transcripts.fa | samtools view -b - > ${sample}_combined.bam

