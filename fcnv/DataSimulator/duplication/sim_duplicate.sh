#!/bin/bash

phase_sites=../../chr20/trio.phase.vcf
bamfile=../../chr20/__M.part.bam
plasmaFile=../../chr20/__plasma.part.bam
readLength=100
plasmaFetusRate=0.13
plasmaCoverage=78

chromosome=chr20
begin=10181440
end=11181440
haplotype=B
source=IM1

region=$chromosome':'$begin'-'$end

pid=$$
logfile=log_simDup.$pid.log
exec > $logfile 2>&1

echo "$region $haplotype $source" 

echo "Copying all the reads in the region..."
samtools view $bamfile $region > raw_reads$pid
# Copying the header
samtools view -H $bamfile $region > inside$pid.sam

echo "Filtering for the correct haplotype..."
python duplication.py raw_reads$pid $phase_sites $source $haplotype > filteredResults$pid

echo "Down sampling the results..."
numReads=`wc -l filteredResults$pid | awk '{print $1}'`
python duplication_down_sampler.py filteredResults$pid $plasmaFetusRate $plasmaCoverage $readLength $numReads $region >> inside$pid.sam
samtools view -S -b inside$pid.sam > inside$pid.bam 

echo "Merging"
# Outside has the headers
samtools merge -f $source-$haplotype-$region-duplicate.bam $plasmaFile inside$pid.bam
echo "Sorting"
samtools sort $source-$haplotype-$region-duplicate.bam $source-$haplotype-$region-duplicate.sort
echo "Indexing"
samtools index $source-$haplotype-$region-duplicate.sort.bam

rm $source-$haplotype-$region-duplicate.bam
rm raw_reads$pid inside$pid.sam filteredResults$pid inside$pid.bam 
