#!/bin/bash

# $1 - maternal reads bam file
# $2 - paternal reads bam file
# $3 - fetal reads bam file
# $4 - reference sequence .fa file
# $5 - region to extract

if [ $# -ne 5 ]; then
    echo "Wrong number of arguments: got $#, expected 5"
    exit
fi

#echo "continued"
#exit
<<comment
# (1) remove PCR duplicates
echo "removing PCR duplicates:"
samtools view -bu $1 $5 | samtools rmdup - - > __M.part.bam &
samtools view -bu $2 $5 | samtools rmdup - - > __P.part.bam &
samtools view -bu $3 $5 | samtools rmdup - - > __F.part.bam &
wait
samtools index __M.part.bam &
samtools index __P.part.bam &
samtools index __F.part.bam &
wait

echo "-------- step 1 done ----------"

# (2) genotype M, P, F

for file in __M.part __P.part __F.part
do
    echo "genotyping $file"
    samtools mpileup -uD -r $5 -f $4 $file.bam | bcftools view -bvcg - > $file.genotype.raw.bcf &
done
wait

for file in __M.part __P.part __F.part
do
    bcftools view $file.genotype.raw.bcf | vcfutils.pl varFilter -D100 > $file.genotype.vcf &
    # ???? what limit for depth of coverage to use?
done
wait

for file in __M.part __P.part __F.part
do
    #extract only SNPs with reasonable quality score
    cat $file.genotype.vcf | ./extract_snps.awk -v qlimit=50 > $file.snps.vcf &
done
wait

echo "-------- step 2 done ----------"
comment

# (3) mix reads to get plasma-like reads
for gnm in M F; do
    echo "getting number of reads and coverage for $gnm"
    count_var=$gnm'_count'
    eval $count_var=$(samtools view __$gnm.part.bam | wc -l)
    
    coverage_var=$gnm'_coverage'
    eval $coverage_var=$(samtools mpileup -D -r $5 __$gnm.part.bam | awk '{sum+=$4; count+=1} END {print sum/count}')
        echo "  >> got ${!count_var} and ${!coverage_var}"
done

target_coverage=80
mixture=0.1

M_multiplier=`echo "scale=5; ($target_coverage * (1.0  - $mixture)) / $M_coverage"|bc`
F_multiplier=`echo "scale=5; ($target_coverage * $mixture) / $F_coverage"|bc`

#for debugging set multipliers by hard
#M_multiplier=0.543
#F_multiplier=0.123
echo "we need $M_multiplier x M and $F_multiplier x F reads" 

#sample reads to simulate plasma reads
temp_file='__temp.sam'
file_count=-1
samtools view -H __M.part.bam > $temp_file
samtools view -H __F.part.bam >> $temp_file

for gnm in M F; do
    frac=$gnm'_multiplier'
    frac=${!frac}
    #echo "+++> $frac"

    while [ $(bc <<< "$frac >= 1") -eq 1 ]; do
        file_count=$(($file_count+1))
        samtools view __$gnm.part.bam > $temp_file$file_count &

        frac=`echo "scale=5; $frac - 1"|bc`
        #echo $frac
    done
    if [ $(bc <<< "$frac > 0") -eq 1 ]; then
        file_count=$(($file_count+1))
        samtools view -s $frac __$gnm.part.bam > $temp_file$file_count &
    fi    
done
wait

#concatenate the temp files
for (( c=0; c<=$file_count; c++ ))
do
    cat $temp_file$c >> $temp_file
done

#sort the sam file: turn it to a bam, sort, convert back to sam
samtools view -Sbu $temp_file | samtools sort - __plasma.sort
rm $temp_file*
plasma_file='plasma.sort.sam'
samtools view -q 20  __plasma.sort.bam > $plasma_file
rm __plasma.sort.bam
echo "-------- step 3 done ----------"

# (4) get nucleotides counts for individual SNP positions (union of the SNP positions found in (2)
# when genotyping the parental genomes)

#plasma_snps='__plasma_mapped.vcf'
#samtools mpileup -uD -f /filer/hg19/hg19.fa plasma.sort.bam | bcftools view -cg - | ./extract_snps.awk -v qlimit=0 > $plasma_snps

echo "combining"
time pypy process_vcfs.py __M.part.snps.vcf __P.part.snps.vcf __F.part.snps.vcf $plasma_file
echo "-------- step 4 done ----------"

# CLEAN-UP
#rm __*.*


