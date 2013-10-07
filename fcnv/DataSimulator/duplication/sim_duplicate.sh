#!/bin/bash

# 1-CNV start position; 2-CNV length; 3-origin(M/P); 4-haplotype(A/B); 5-path to data;
# 6-path where to store the simulated plasma BAM; 7-path to dir with scripts

data_path=$5
results_path=$6
exec_path=$7
phase_sites=$data_path/trio.phase.vcf
bamfile=$data_path'/__'$3'.part.bam'
plasmaFile=$data_path/__plasma.part.bam
tmp_path=/tmp
readLength=100
plasmaFetusRate=0.13

chromosome=chr20
begin=$1
end=$(($1 + $2))
haplotype=$4
source='I'$3'1'

region=$chromosome':'$begin'-'$end

pid=$$
logfile=$results_path/log_simDup.$pid.log
exec > $logfile 2>&1


echo "SimDuplication: $region $haplotype $source $3"
date

#temp files names
tmp_pileup_file=$tmp_path/plasma_pileup$pid.txt
raw_reads_file=$tmp_path/raw_reads$pid.sam
inside_sam_file=$tmp_path/inside$pid.sam
inside_bam_file=$tmp_path/inside$pid.bam
filtered_res_file=$tmp_path/filteredResults$pid.sam
plasma_file_prefix=$tmp_path/$source-$haplotype-$region-duplicate

echo "Pileuping plasma reads in the region..."
tmp_region=$chromosome':'$(($begin - 1000))'-'$(($end + 1000))
samtools mpileup $plasmaFile -q10 -Q10 -r $tmp_region | awk '{print $2, $4}' > $tmp_pileup_file

echo "Copying all the reads in the region..."
samtools view $bamfile $region > $raw_reads_file
# Copying the header
samtools view -H $bamfile $region > $inside_sam_file

echo "Filtering for the correct haplotype..."
pypy $exec_path/duplication.py $raw_reads_file $phase_sites $haplotype > $filtered_res_file

echo "Down sampling the results..."
numReads=`wc -l $filtered_res_file | awk '{print $1}'`
pypy $exec_path/duplication_down_sampler.py $filtered_res_file $numReads $readLength $tmp_pileup_file $plasmaFetusRate $region >> $inside_sam_file
samtools view -S -b $inside_sam_file > $inside_bam_file

echo "Merging"
# Outside has the headers
samtools merge -f $plasma_file_prefix.bam $plasmaFile $inside_bam_file
echo "Sorting"
samtools sort $plasma_file_prefix.bam $plasma_file_prefix.sort
echo "Indexing"
samtools index $plasma_file_prefix.sort.bam

#move to results_path
mv $plasma_file_prefix.sort.bam* $results_path/

rm $plasma_file_prefix.bam
rm $tmp_pileup_file $raw_reads_file $filtered_res_file $inside_sam_file $inside_bam_file 

echo "sim_duplicate script - DONE."

