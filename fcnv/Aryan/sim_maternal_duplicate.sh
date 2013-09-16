bamfile=~/fetusProject/I1/chr20/I1_M_chr20.bam
region=chr20:10181440-10281440
#phase_sites=~/fetusProject/I1/chr20/i1_phase_blocks_and_sites.txt 
phase_sites=./trio.phase.vcf
haplotype=B

echo "copying all the reads in the region..."
samtools view $bamfile $region > tmp_reads
echo "done."

echo "filtering for the correct haplotype"
python maternal_duplicate.py tmp_reads $phase_sites $haplotype > resultsB 
echo "done."

