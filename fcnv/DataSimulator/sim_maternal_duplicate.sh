bamfile=~/fetusProject/I1/chr20/I1_M_chr20.bam
region=chr20:10181440-10281440
#phase_sites=~/fetusProject/I1/chr20/i1_phase_blocks_and_sites.txt 
phase_sites=./trio.phase.vcf
haplotype=B
source=IM1
plasmaFetusRate=0.13
plasmaCoverage=78
readLength=100

echo "Copying all the reads in the region..."
samtools view $bamfile $region > raw_reads
echo "done."

echo "Filtering for the correct haplotype..."
python maternal_duplicate.py raw_reads $phase_sites $source $haplotype > filteredResults
echo "done."

echo "Down sampling the results..."
numReads=`wc -l filteredResults | awk '{print $1}'`
python maternal_duplicate_down_sampler.py filteredResults $plasmaFetusRate $plasmaCoverage $readLength $numReads $region > $source-$haplotype-$region-duplicate.sam
echo "done."


