#!/bin/bash

plasmaFetusRate=$1
chrom=$2
phase_sites=/dupa-filer/laci/I1/$chrom/trio.phase.vcf
bamfile=/dupa-filer/laci/I1/$chrom/__plasma.part.bam
samfile=/tmp/__plasma.$chrom.sam
output_bamfile=/tmp/'__plasma-'$chrom'-'$plasmaFetusRate'rate.sort.bam'
output_samfile=/tmp/'__plasma-'$chrom'-'$plasmaFetusRate'rate.sort.sam'

# Copying the header
samtools view -H $bamfile -o $output_samfile
samtools view $bamfile -l /dupa-filer/laci/I1/$chrom/regions.bed -o $samfile

echo "Filtering the reads to reduce fetus DNA rate"
pypy fetusRateAdjuster.py $samfile $phase_sites $plasmaFetusRate >> $output_samfile
samtools view -S -b $output_samfile -o $output_bamfile

echo "Indexing"
samtools index $output_bamfile
mv $output_bamfile* /dupa-filer/laci/I1/$chrom/

echo "sim_fetusRate script - DONE."

