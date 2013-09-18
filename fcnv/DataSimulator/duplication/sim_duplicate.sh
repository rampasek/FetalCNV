phase_sites=./trio.phase.vcf
bamfile=~/fetusProject/I1/chr20/I1_M_chr20.bam
plasmaFile=
readLength=100
plasmaFetusRate=0.13
plasmaCoverage=78

chromosome=chr20
begin=10181440
end=10281440
haplotype=B
source=IM1

region=$chromosome':'$begin'-'$end

echo "Copying all the reads in the region..."
samtools view $bamfile $region > raw_reads
# Copying the header
samtools view -H $bamfile $region > inside.sam

echo "Filtering for the correct haplotype..."
python duplication.py raw_reads $phase_sites $source $haplotype > filteredResults

echo "Down sampling the results..."
numReads=`wc -l filteredResults | awk '{print $1}'`
python duplication_down_sampler.py filteredResults $plasmaFetusRate $plasmaCoverage $readLength $numReads $region >> inside.sam
samtools view -S -b inside.sam > inside.bam 

echo "Merging"
# Outside has the headers
samtools merge -f $source-$haplotype-$region-duplicate.bam $plasmaFile inside.bam
samtools sort $source-$haplotype-$region-duplicate.bam $source-$haplotype-$region-duplicate.sort
samtools index $source-$haplotype-$region-duplicate.sort.bam

rm $source-$haplotype-$region-duplicate.bam
