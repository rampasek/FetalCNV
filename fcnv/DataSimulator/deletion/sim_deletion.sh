bamfile=~/fetusProject/I1/chr20/I1_M_chr20.bam
region=chr20:10181440-10281440
regionCompliment='chr20:1-10181339 chr20:10281541'
#phase_sites=~/fetusProject/I1/chr20/i1_phase_blocks_and_sites.txt 
phase_sites=./trio.phase.vcf
haplotype=B
source=IM1
plasmaFetusRate=0.13
plasmaCoverage=78
readLength=100

echo "Copying all the reads in the region..."
samtools view $bamfile $region > raw_reads
samtools view -H $bamfile $region > $haplotype-$region-delete.sam

echo "Filtering the reads that are not deleted..."
python deletion.py raw_reads $phase_sites $plasmaFetusRate $haplotype >> $haplotype-$region-delete.sam
samtools view -S -b $haplotype-$region-delete.sam > inside.bam 

echo "Adding reads located out of the region..."
samtools view -h -b $bamfile $regionCompliment > outside.bam

echo "Merging"
# Outside has the headers
samtools merge $haplotype-$region-delete.bam outside.bam inside.bam

