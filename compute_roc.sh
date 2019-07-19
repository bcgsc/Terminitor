#!/bin/bash

prefix=$1
touch $2

for i in $(seq 1 1 1)
  do 
    awk -v thred=$i '$5>=thred' ${prefix}.forward.cluster.bed > t1
    awk -v thred=$i '$5>=thred' ${prefix}.reverse.cluster.bed > t2
    # echo $i
    # /projects/cheny_prj/KLEAT_benchmarking/scripts/compute_metrics.py --i1 t1 --i2 t2 --r1 /projects/cheny_prj/KLEAT_benchmarking/polyA_seq/UHR.combined_noChrM.forward.bed --r2 /projects/cheny_prj/KLEAT_benchmarking/polyA_seq/UHR.combined_noChrM.reverse.bed >> $2
    /projects/btl_scratch/cheny/NN_paper/pipeline/compute_metrics.py --i1 t1 --i2 t2 --r1 $3 --r2 $4 $5 >> $2
    # /projects/cheny_prj/KLEAT_benchmarking/scripts/compute_metrics.py --i1 t1 --i2 t2 --r1 ../../../Ensembl_havana_with_polyAtail.forward.cluster.bed --r2 ../../../Ensembl_havana_with_polyAtail.reverse.cluster.bed >> $2
  done


# for i in $(seq 0.991 0.001 0.999)
  # do
    #awk -v thred=$i '$5>=thred' ${prefix}.forward.cluster.bed > t1
    #awk -v thred=$i '$5>=thred' ${prefix}.reverse.cluster.bed > t2
    # echo $i
    # /projects/cheny_prj/KLEAT_benchmarking/scripts/compute_metrics.py --i1 t1 --i2 t2 --r1 /projects/cheny_prj/KLEAT_benchmarking/polyA_seq/UHR.combined_noChrM.forward.bed --r2 /projects/cheny_prj/KLEAT_benchmarking/polyA_seq/UHR.combined_noChrM.reverse.bed >> $2
    #/projects/btl_scratch/cheny/NN_paper/pipeline/compute_metrics.py --i1 t1 --i2 t2 --r1 $3 --r2 $4 $5 >> $2
    # /projects/cheny_prj/KLEAT_benchmarking/scripts/compute_metrics.py --i1 t1 --i2 t2 --r1 ../../../Ensembl_havana_with_polyAtail.forward.cluster.bed --r2 ../../../Ensembl_havana_with_polyAtail.reverse.cluster.bed >> $2
  #done
