#!/bin/bash

#SBATCH --ntasks=48

script_dir='/projects/btl_scratch/cheny/NN_paper/model/'
train_dir=$4

${script_dir}/three_label_classification.py ${train_dir}/PolyA_DB3_positive.fasta ${train_dir}/Ensembl_cs_trans_with_no_polyAtail.fasta ${train_dir}/Four_db_negative.fasta $1 $2 $3

# ${script_dir}/three_label_classification.py /projects/btl_scratch/cheny/NN/length_experiment/training_set/100_100/PolyA_DB3_positive_40000.fa /projects/btl_scratch/cheny/NN/length_experiment/training_set/100_100/Ensembl_cs_trans_with_no_polyAtail_40000.fa /projects/btl_scratch/cheny/NN/length_experiment/training_set/100_100/Four_db_negative_40000_6.fa $1 $2 $3 $4 $5
