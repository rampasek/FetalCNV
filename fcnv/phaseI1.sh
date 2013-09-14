#!/bin/bash

# $1 - maternal reads bam file
# $2 - paternal reads bam file
# $3 - fetal reads bam file
# $4 - plasma reads bam file
# $5 - reference sequence .fa file
# $6 - region to extract

logfile=log_$6.$$.log
exec > $logfile 2>&1

date

if [ $# -ne 6 ]; then
    echo "Wrong number of arguments: got $#, expected 6"
    exit
fi

bindir="/dupa-filer/laci/bin"

#parse out chromosome ID
region=$6
chr=`echo $region | awk 'BEGIN {FS=":"}{print $1}'`
reference=$5

<<comment
# (1) remove PCR duplicates
echo "merging and removing PCR duplicates:"
samtools merge -R $region - $1 | samtools rmdup - __M.part.bam &
samtools merge -R $region - $2 | samtools rmdup - __P.part.bam &
samtools merge -R $region - $3 | samtools rmdup - __F.part.bam &
samtools merge -R $region - $4 | samtools rmdup - __plasma.part.bam &
#samtools view -bu $4 $region | samtools rmdup - - > __plasma.part.bam &
wait
samtools index __M.part.bam &
samtools index __P.part.bam &
samtools index __F.part.bam &
samtools index __plasma.part.bam &
wait
samtools view -h __M.part.bam > __M.part.sam &
samtools view -h __P.part.bam > __P.part.sam &
#wait
comment
echo "-------- step 1 done ----------"

# (2) genotype M, P, F, filter and phase

for file in __M.part __P.part __F.part
do
    echo "genotyping $file"
    samtools mpileup -uD -r $region -f $reference $file.bam | bcftools view -bvcg - > $file.genotype.raw.bcf &
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
    #annotate SNPs by rs-ids from dbSNP
    echo  "annotating $file"
    java -jar ~/apps/snpEff/SnpSift.jar annotate -v /dupa-filer/laci/dbSnp.vcf $file.genotype.vcf > $file.genotype.annot.vcf &
done
wait

for file in __M.part __P.part __F.part
do
    #extract only SNPs with reasonable quality score
    cat $file.genotype.annot.vcf | extract_annot_snps.awk -v qlimit=50 > $file.snps.annot.vcf &
done
wait

#----------------
# filter out SNP positions for which we don't have confident calls (or homoz. ref.) evidance in both M and P
echo "Filtering SNP positions"
time $bindir/filter_vcfs.py __M.part.snps.annot.vcf __P.part.snps.annot.vcf __M.part.sam __P.part.sam


for file in __M.part __P.part 
do
    #bgzip and tabix index the VCFs
    bgzip $file.snps.annot.ftr.vcf
    tabix -p vcf $file.snps.annot.ftr.vcf.gz
done
bgzip __F.part.snps.annot.vcf
tabix -p vcf __F.part.snps.annot.vcf.gz

#merge VCFs
vcf-merge __M.part.snps.annot.vcf.gz __P.part.snps.annot.vcf.gz __F.part.snps.annot.vcf.gz

refpanels=/dupa-filer/laci/reference_panels/$chr.1kg.ref.phase1_release_v3.20101123.vcf.gz
for prefix in __M __P
do
    #convert annotated VCF to BEAGLE format (creates $file.bgl.gz, .markers, and .int)
    cat $prefix.part.snps.annot.ftr.vcf | java -jar ~/apps/jar/vcf2beagle.jar ? $prefix
    #unify alleles for the markers in $file.markers
    #pypy $bindir/union_markers.py $prefix.markers $refpanels.markers > $prefix.mod.markers
    #phase haplotypes
    java -Xmx10000m -jar ~/apps/jar/beagle.jar markers=$prefix.mod.markers phased=$refpanels  unphased=$prefix.bgl.gz missing=? out=$prefix omitprefix=true &
done
wait
        
for prefix in __M __P
do
    #convert haplotypes from BEAGLE to VCF format
    java -jar ~/apps/jar/beagle2vcf.jar $chr $prefix.markers $prefix.bgl.gz.phased.gz ? > $prefix.phased.vcf &
done
wait


# CLEAN-UP
#rm __*.*



