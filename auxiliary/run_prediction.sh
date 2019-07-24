#!/bin/bash
#SBATCH --ntasks=48
#SBATCH -N 1-1
#SBATCH --exclude='n[121-124,129-132,224-228]'
#SBATCH --mem=360G

/usr/bin/time -v ./test_model.py ../extracted_reads/${1}.fa ../../length_experiments/mouse/model/100_100/pred_model1_best_weights.hdf5 200 > $1

/usr/bin/time -v ./test_model.py ../extracted_reads/${2}.fa ../../length_experiments/mouse/model/100_100/pred_model1_best_weights.hdf5 200 > $2

