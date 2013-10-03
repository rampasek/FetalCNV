#!/bin/bash

ref=$1
pos=$2
reads1=$3
reads2=$4
reads3=$5
prefix=$6
valid_regions=$7
pile_prefix=$8

for num in 1 2 3
do
    reads_file_var=reads$num
    reads_file=${!reads_file_var}
    #echo "samtools mpileup -f $ref -uDSI -Q10 -q10 -l $pos $reads_file | bcftools view -g - > $prefix.$num.vcf"
    samtools mpileup -f $ref -uDSI -Q20 -q10 -l $pos $reads_file | bcftools view -g - > $prefix.$num.vcf &
done

#pileup plasma and maternal reads to get DOC for all positions
for num in 1 2
do
    reads_file_var=reads$num
    reads_file=${!reads_file_var}
    samtools mpileup -Q20 -q10 -l $valid_regions $reads_file | awk '{print $2,$4}' > $pile_prefix.$num.txt &
done

wait
