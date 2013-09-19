#!/bin/bash

phase_sites=../../chr20/trio.phase.vcf
bamfile=../../chr20/__plasma.part.bam
readLength=100
plasmaFetusRate=0.13

chromosome=chr20
begin=$1
end=$(($1 + $2))
haplotype=B

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
python deletion.py raw_reads$pid $phase_sites $plasmaFetusRate $haplotype >> inside$pid.sam
samtools view -S -b inside$pid.sam > inside$pid.bam 

echo "Adding reads located out of the region..."
samtools view -h -b $bamfile $regionCompliment > outside$pid.bam

echo "Merging"
# Outside has the headers
samtools merge -f $haplotype-$region-delete.bam outside$pid.bam inside$pid.bam
echo "Sorting"
samtools sort $haplotype-$region-delete.bam $haplotype-$region-delete.sort
echo "Indexing"
samtools index $haplotype-$region-delete.sort.bam

rm $haplotype-$region-delete.bam
rm outside$pid.bam inside$pid.sam inside$pid.bam raw_reads$pid

