#!/bin/bash

ref=$1
pos=$2
reads1=$3
reads2=$4
reads3=$5
tmp_prefix=$6
wbegin=$7
wend=$8
chr=$9
pile_prefix=${10}

for num in 1 2 3
do
    reads_file_var=reads$num
    reads_file=${!reads_file_var}
    #echo "samtools mpileup -f $ref -uDSI -Q10 -q10 -l $pos $reads_file | bcftools view -g - > $prefix.$num.vcf"
    samtools mpileup -f $ref -uDSI -Q10 -q10 -l $pos $reads_file | bcftools view -g - > $tmp_prefix.$num.vcf &
done

nums='1' #'1 2'
#pileup plasma and maternal reads to get DOC for all positions
for num in $nums
do
    echo "pileup plasma: $chr:$wbegin-$wend"
    reads_file_var=reads$num
    reads_file=${!reads_file_var}
    samtools mpileup -Q10 -q10 -r $chr:$wbegin-$wend $reads_file | awk '{print $2,$4}' > $pile_prefix.$num.txt & #$tmp_prefix.$num.txt
done
wait


#move the pileup files to their destination
#for num in $nums
#do
#    mv $tmp_prefix.$num.txt $pile_prefix.$num.txt
#done

