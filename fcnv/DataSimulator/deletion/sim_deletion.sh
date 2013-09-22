#!/bin/bash

# 1-CNV start position; 2-CNV length; 3-haplotype(A/B); 4-path to data;
# 5-path where to store the simulated plasma BAM; 6-path to dir with scripts

data_path=$4
results_path=$5
exec_path=$6
phase_sites=$data_path/trio.phase.vcf
bamfile=$data_path/__plasma.part.bam
readLength=100
plasmaFetusRate=0.13

chromosome=chr20
begin=$1
end=$(($1 + $2))
haplotype=$3

regionCompliment=$chromosome':1-'$((begin-readLength))' '$chromosome':'$((end+readLength))
region=$chromosome':'$begin'-'$end

pid=$$
logfile=log_simDel.$pid.log
exec > $logfile 2>&1

echo "$region $haplotype $source"

echo "Copying all the reads in the region..."
samtools view $bamfile $region > raw_reads$pid
# Copying the header
samtools view -H $bamfile $region > inside$pid.sam

echo "Filtering the reads that are not deleted..."
python $exec_path/deletion.py raw_reads$pid $phase_sites $plasmaFetusRate $haplotype >> inside$pid.sam
samtools view -S -b inside$pid.sam > inside$pid.bam

echo "Adding reads located out of the region..."
samtools view -h -b $bamfile $regionCompliment > outside$pid.bam

echo "Merging"
# Outside has the headers
samtools merge -f $haplotype-$region-delete.bam outside$pid.bam inside$pid.bam
echo "Sorting"
samtools sort $haplotype-$region-delete.bam $haplotype-$region-delete.sort
echo "Indexing"
samtools index $results_path/$haplotype-$region-delete.sort.bam

#move to results_path
mv $source-$haplotype-$region-duplicate.sort.bam* $results_path/

rm $haplotype-$region-delete.bam
rm outside$pid.bam inside$pid.sam inside$pid.bam raw_reads$pid

echo "sim_deletion script - DONE."

