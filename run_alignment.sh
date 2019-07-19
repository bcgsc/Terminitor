#!/bin/bash
#SBATCH --ntasks=12

/usr/bin/time -v minimap2 -ax splice hg38.mmi -t 12 ../assembly/UHRR-C1_S1/rnabloom.transcripts.fa | samtools view -b - > UHRR-C1-S1.bam

/usr/bin/time -v minimap2 -ax splice hg38.mmi -t 12 ../assembly/UHRR-C2_S2/rnabloom.transcripts.fa | samtools view -b - > UHRR-C2-S2.bam

/usr/bin/time -v minimap2 -ax splice hg38.mmi -t 12 ../assembly/HBRR-C4_S3/rnabloom.transcripts.fa | samtools view -b - > HBRR-C4-S3.bam

/usr/bin/time -v minimap2 -ax splice hg38.mmi -t 12 ../assembly/HBRR-C6_S4/rnabloom.transcripts.fa | samtools view -b - > HBRR-C6-S4.bam


