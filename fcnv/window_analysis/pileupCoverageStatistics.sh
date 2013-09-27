bamfile=../Data/__plasma.part.bam
bedfile=regions.bed
samtools mpileup $bamfile -l $bedfile | awk '{print $2,$4}'
