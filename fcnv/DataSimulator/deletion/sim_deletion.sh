#!/bin/bash

# 1-CNV start position; 2-CNV length; 3-haplotype(A/B); 4-path to data;
# 5-path where to store the simulated plasma BAM; 6-path to dir with scripts

data_path=$4
results_path=$5
exec_path=$6
phase_sites=$data_path/trio.phase.vcf
bamfile=$data_path/__plasma.part.bam
tmp_path=/tmp
readLength=100
plasmaFetusRate=0.13

chromosome=chr20
begin=$1
end=$(($1 + $2))
haplotype=$3

regionCompliment=$chromosome':1-'$((begin-readLength))' '$chromosome':'$((end+readLength))
region=$chromosome':'$begin'-'$end

pid=$$
logfile=$results_path/log_simDel.$pid.log
exec > $logfile 2>&1

echo "SimDeletion: $region $haplotype $source"

#temp files names
raw_reads_file=$tmp_path/raw_reads$pid
inside_sam_file=$tmp_path/inside$pid.sam
inside_bam_file=$tmp_path/inside$pid.bam
outside_bam_file=$tmp_path/outside$pid.bam
plasma_file_prefix=$tmp_path/$haplotype-$region-delete

echo "Copying all the reads in the region..."
samtools view $bamfile $region > $raw_reads_file
# Copying the header
samtools view -H $bamfile $region > $inside_sam_file

echo "Filtering the reads that are not deleted..."
pypy $exec_path/deletion.py $raw_reads_file $phase_sites $plasmaFetusRate $haplotype >> $inside_sam_file
samtools view -S -b $inside_sam_file > $inside_bam_file

echo "Adding reads located out of the region..."
samtools view -h -b $bamfile $regionCompliment > $outside_bam_file

echo "Merging"
# Outside has the headers
samtools merge -f $plasma_file_prefix.bam $outside_bam_file $inside_bam_file
echo "Sorting"
samtools sort $plasma_file_prefix.bam $plasma_file_prefix.sort
echo "Indexing"
samtools index $plasma_file_prefix.sort.bam

#move to results_path
mv $plasma_file_prefix.sort.bam* $results_path/

rm $plasma_file_prefix.bam
rm $outside_bam_file $inside_sam_file $inside_bam_file $raw_reads_file

echo "sim_deletion script - DONE."

