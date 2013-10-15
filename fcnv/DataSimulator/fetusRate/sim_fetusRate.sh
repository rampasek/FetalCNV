#!/bin/bash

plasmaFetusRate=$1
phase_sites=trio.phase.vcf
bamfile=__plasma.part.bam
samfile=__plasma.sam
output_bamfile='__plasma-'$plasmaFetusRate'rate.sort.bam'
output_samfile='__plasma-'$plasmaFetusRate'rate.sort.sam'

# Copying the header
samtools view -H $bamfile -o $output_samfile

echo "Filtering the reads to reduce fetus DNA rate"
pypy fetusRateAdjuster.py $samfile $phase_sites $plasmaFetusRate >> $output_samfile
samtools view -S -b $output_samfile -o $output_bamfile

echo "Indexing"
samtools index $output_bamfile

echo "sim_fetusRate script - DONE."

