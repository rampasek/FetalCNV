phase_sites=./trio.phase.vcf
bamfile=~/fetusProject/I1/chr20/I1_M_chr20.bam
readLength=100
plasmaFetusRate=0.13

chromosome=chr20
begin=10181440
end=10281440
haplotype=B

regionCompliment=$chromosome':1-'$((begin-readLength))' '$chromosome':'$((end+readLength))
region=$chromosome':'$begin'-'$end

echo "Copying all the reads in the region..."
samtools view $bamfile $region > raw_reads
# Copying the header
samtools view -H $bamfile $region > inside.sam

echo "Filtering the reads that are not deleted..."
python deletion.py raw_reads $phase_sites $plasmaFetusRate $haplotype >> inside.sam
samtools view -S -b inside.sam > inside.bam 

echo "Adding reads located out of the region..."
samtools view -h -b $bamfile $regionCompliment > outside.bam

echo "Merging"
# Outside has the headers
samtools merge -f $haplotype-$region-delete.bam outside.bam inside.bam

