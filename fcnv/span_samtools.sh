#!/bin/bash

ref=$1
pos=$2
reads1=$3
reads2=$4
reads3=$5
prefix=$6

for num in 1 2 3
do
    reads_file=reads$num
    echo "samtools mpileup -f $ref -u -Q15 -q10 -l $pos ${!reads_file} | bcftools view -g - > $prefix.$num.vcf"
    samtools mpileup -f $ref -uDSI -Q15 -q10 -l $pos ${!reads_file} | bcftools view -g - > $prefix.$num.vcf &
done
wait
