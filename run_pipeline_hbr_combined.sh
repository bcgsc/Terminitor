#!/bin/bash
#SBATCH --ntasks=48
#SBATCH -N 1-1
#SBATCH --exclude='n[121-124,129-132,224-228]'
#SBATCH --mem=360G

/usr/bin/time -v /home/kmnip/jdk1.8.0_101/jre/bin/java -jar /projects/kmnip_prj/RNA-Bloom/dist/RNA-Bloom_v1.0.0/RNA-Bloom.jar -left /projects/btl2/kmnip/rna-bloom/example/hbr/mRNA-Brain-C4_S3_L001_R2_001.fastq.gz /projects/btl2/kmnip/rna-bloom/example/hbr/mRNA-Brain-C4_S3_L002_R2_001.fastq.gz /projects/btl2/kmnip/rna-bloom/example/hbr/mRNA-Brain-C6_S4_L001_R2_001.fastq.gz /projects/btl2/kmnip/rna-bloom/example/hbr/mRNA-Brain-C6_S4_L002_R2_001.fastq.gz -right /projects/btl2/kmnip/rna-bloom/example/hbr/mRNA-Brain-C4_S3_L001_R1_001.fastq.gz /projects/btl2/kmnip/rna-bloom/example/hbr/mRNA-Brain-C4_S3_L002_R1_001.fastq.gz /projects/btl2/kmnip/rna-bloom/example/hbr/mRNA-Brain-C6_S4_L001_R1_001.fastq.gz /projects/btl2/kmnip/rna-bloom/example/hbr/mRNA-Brain-C6_S4_L002_R1_001.fastq.gz -revcomp-right -t 12 -outdir HBRR_combined -a 4 -e 1 -stratum 01 -ss -ntcard -fpr 0.005

