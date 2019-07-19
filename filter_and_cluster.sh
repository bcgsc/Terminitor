#!/bin/bash

pred=$1
fa=$2
prefix=$3
curdir=`pwd`

paste <(grep '>' $fa) $pred | cut -f 2 -d '>' > ${prefix}.prediction

awk -v FS='_' -v OFS='\t' '{print $2,$3,$3+1,$1,$4}' ${prefix}.prediction | awk -v OFS='\t' '{print $1,$2,$3,$4,$6,$5}' > ${prefix}.bed
sed -i 's/R/-/g' ${prefix}.bed
sed -i 's/F/+/g' ${prefix}.bed

/gsc/btl/linuxbrew/bin/sort -k1,1 -k2,2n -k6,6 ${prefix}.bed > ${prefix}.sorted.bed

awk '$6=="+"' ${prefix}.sorted.bed > ${prefix}.forward.bed
awk '$6=="-"' ${prefix}.sorted.bed > ${prefix}.reverse.bed
/projects/btl_scratch/cheny/NN/scripts/cluster.py -d 31 -o $prefix.forward ${prefix}.forward.bed
/projects/btl_scratch/cheny/NN/scripts/cluster.py -d 31 -o $prefix.reverse ${prefix}.reverse.bed

rm ${prefix}.bed ${prefix}.prediction 
rm ${prefix}.forward.bed ${prefix}.reverse.bed
