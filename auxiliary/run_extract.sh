#!/bin/bash
#SBATCH --ntasks=48
#SBATCH -N 1-1
#SBATCH --exclude='n[121-124,129-132,224-228]'
#SBATCH --mem=360G

sample=$1
/usr/bin/time -v ./extract_from_sequences.py ../../Annotation/Homo_sapiens.GRCh38.94.transcripts.gtf ../../Annotation/Homo_sapiens.GRCh38.94.chr.gtf ../alignment/${sample}_combined.bam ../../Annotation/Homo_sapiens.GRCh38.dna.primary_assembly.fa ${sample}_combined_uniq.fa 100 100

