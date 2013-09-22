#!/bin/bash

# 1-CNV start position; 2-CNV length; 3-origin(M/P); 4-haplotype(A/B); 5-path to data;
# 6-path where to store the simulated plasma BAM; 7-path to dir with scripts

data_path=$5
results_path=$6
exec_path=$7
phase_sites=$data_path/trio.phase.vcf
bamfile=$data_path'/__'$3'.part.bam'
plasmaFile=$data_path/__plasma.part.bam
readLength=100
plasmaFetusRate=0.13
plasmaCoverage=78

chromosome=chr20
begin=$1
end=$(($1 + $2))
haplotype=$4
source='I'$3'1'

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
python $exec_path/duplication.py raw_reads$pid $phase_sites $source $haplotype > filteredResults$pid

echo "Down sampling the results..."
numReads=`wc -l filteredResults$pid | awk '{print $1}'`
python $exec_path/duplication_down_sampler.py filteredResults$pid $plasmaFetusRate $plasmaCoverage $readLength $numReads $region >> inside$pid.sam
samtools view -S -b inside$pid.sam > inside$pid.bam

echo "Merging"
# Outside has the headers
samtools merge -f $source-$haplotype-$region-duplicate.bam $plasmaFile inside$pid.bam
echo "Sorting"
samtools sort $source-$haplotype-$region-duplicate.bam $source-$haplotype-$region-duplicate.sort
echo "Indexing"
samtools index $results_path/$source-$haplotype-$region-duplicate.sort.bam

#move to results_path
mv $source-$haplotype-$region-duplicate.sort.bam* $results_path/

rm $source-$haplotype-$region-duplicate.bam
rm raw_reads$pid inside$pid.sam filteredResults$pid inside$pid.bam

echo "sim_duplicate script - DONE."
