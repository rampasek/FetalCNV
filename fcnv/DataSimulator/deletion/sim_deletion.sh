#!/bin/bash

phase_sites=../../chr20/trio.phase.vcf
bamfile=../../chr20/__plasma.part.bam
readLength=100
plasmaFetusRate=0.13

chromosome=chr20
begin=10181440
end=20181440
haplotype=B

regionCompliment=$chromosome':1-'$((begin-readLength))' '$chromosome':'$((end+readLength))
region=$chromosome':'$begin'-'$end

logfile=log_simDel.$$.log
exec > $logfile 2>&1

echo "$region $haplotype $source"

echo "Copying all the reads in the region..."
samtools view $bamfile $region > raw_reads$$    
# Copying the header
samtools view -H $bamfile $region > inside$$.sam

echo "Filtering the reads that are not deleted..."
python deletion.py raw_reads$$ $phase_sites $plasmaFetusRate $haplotype >> inside$$.sam
samtools view -S -b inside$$.sam > inside$$.bam 

echo "Adding reads located out of the region..."
samtools view -h -b $bamfile $regionCompliment > outside$$.bam

echo "Merging"
# Outside has the headers
samtools merge -f $haplotype-$region-delete.bam outside$$.bam inside$$.bam
echo "Sorting"
samtools sort $haplotype-$region-delete.bam $haplotype-$region-delete.bam.sort
echo "Indexing"
samtools index $haplotype-$region-delete.bam.sort.bam

rm $haplotype-$region-delete.bam
rm outside$$.bam inside$$.bam raw_reads$$

