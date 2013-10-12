#!/bin/bash

# 1-chromosome; 2-CNV start position; 3-CNV length; 4-haplotype(A/B); 5-path to data;
# 6-path where to store the simulated plasma BAM; 7-path to dir with scripts

data_path=$5
results_path=$6
exec_path=$7
phase_sites=$data_path/trio.phase.vcf
bamfile=$data_path/__plasma.part.bam
local_bam=/tmp/data/__plasma.part.bam
readLength=100
plasmaFetusRate=0.13
tmp_path=/tmp
if [[ -n "$TMPDIR" ]]; then
    tmp_path=$TMPDIR
fi
if [[ -r "$local_bam" ]]; then
    bamfile=$local_bam
fi

chromosome=$1
begin=$2
end=$(($2 + $3))
haplotype=$4

pid=$$
logfile=$results_path/log_simDel.$pid.log
exec > $logfile 2>&1

regionComplementA=$chromosome':1-'$((begin-readLength))
regionComplementB=$chromosome':'$((end+readLength))
region=$chromosome':'$begin'-'$end

echo "SimDeletion: $region $haplotype $source"
date

#temp files names
raw_reads_file=$tmp_path/raw_reads$pid.sam
inside_sam_file=$tmp_path/inside$pid.sam
inside_bam_file=$tmp_path/inside$pid.bam
outsideA_bam_file=$tmp_path/outsideA$pid.bam
outsideB_bam_file=$tmp_path/outsideB$pid.bam
plasma_file_prefix=$tmp_path/$haplotype-$region-delete

echo "Copying all the reads in the region..."
samtools view $bamfile $region -o $raw_reads_file
# Copying the header
samtools view -H $bamfile $region -o $inside_sam_file

echo "Filtering the reads that are not deleted..."
pypy $exec_path/deletion.py $raw_reads_file $phase_sites $plasmaFetusRate $haplotype >> $inside_sam_file
samtools view -S -b $inside_sam_file -o $inside_bam_file

echo "Adding reads located out of the region..."
samtools view -h -b $bamfile $regionComplementA -o $outsideA_bam_file
samtools view -h -b $bamfile $regionComplementB -o $outsideB_bam_file

echo "Merging"
# Outside has the headers
samtools merge -f $plasma_file_prefix.sort.bam $outsideA_bam_file $inside_bam_file $outsideB_bam_file
#echo "Sorting"
#samtools sort $plasma_file_prefix.bam $plasma_file_prefix.sort
echo "Indexing"
samtools index $plasma_file_prefix.sort.bam

#move to results_path
mv $plasma_file_prefix.sort.bam* $results_path/

#rm $plasma_file_prefix.bam
rm $outsideA_bam_file $outsideB_bam_file $inside_sam_file $inside_bam_file $raw_reads_file

echo "sim_deletion script - DONE."

