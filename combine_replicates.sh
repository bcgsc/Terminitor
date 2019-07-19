#!/bin/bash

l1=$1
l2=$2
prefix=$3

cat $l1 $l2 | /gsc/btl/linuxbrew/bin/sort -k1,1 -k2,2n -k6,6 > ${prefix}.sorted.bed

# /projects/btl_scratch/cheny/NN/scripts/cluster.py -d 31 -o $prefix ${prefix}.sorted.bed
awk '$6=="+"' ${prefix}.sorted.bed > ${prefix}.forward.bed
awk '$6=="-"' ${prefix}.sorted.bed > ${prefix}.reverse.bed
/projects/btl_scratch/cheny/NN/scripts/cluster.py -d 31 -o $prefix.forward ${prefix}.forward.bed
/projects/btl_scratch/cheny/NN/scripts/cluster.py -d 31 -o $prefix.reverse ${prefix}.reverse.bed

rm ${prefix}.forward.bed ${prefix}.reverse.bed

