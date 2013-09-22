#!/bin/bash

data_path='/dupa-filer/laci/I1/chr20/'
plasma_path='/dupa-filer/laci/I1/chr20/sim_plasma/'
results_path='/dupa-filer/laci/I1/chr20/fcnv_data/'
exec_path='/dupa-filer/laci/bin/'

echo "Starting CNV plasma files simulations"
#time ./sim_duplicate.sh 10000000 100000 P &

length_set='10000000 1000000 500000 100000'
haplo_set='A B'
src_set='P M'
CHR_LEN=100000000

for i in {1..3}
do
	for length in $length_set
	do
		for haplo in $haplo_set
		do
			begin=$RANDOM
			begin=$(($begin*$RANDOM))
			begin=$(($begin%$CHR_LEN))

			#qsub -q all.q -R y -pe parallel 6 -l h_vmem=10G -l h_rt=05:00:00 -S $exec_path/sim_deletion.sh $begin $length $haplo $data_path $plasma_path $exec_path
			echo $begin' '$length' '$haplo

			for src in $src_set
			do
				#qsub -q all.q -R y -pe parallel 6 -l h_vmem=10G -l h_rt=05:00:00 -S $exec_path/sim_duplicate.sh $begin $length $src $haplo $data_path $plasma_path $exec_path
				echo $begin' '$length' '$src' '$haplo
			done
		done
	done
done

exit

$exec_path/sim_duplicate.sh 10000000 100000 P B $data_path $plasma_path $exec_path &
$exec_path/sim_deletion.sh 10000000 100000 A $data_path $plasma_path $exec_path &
wait
echo "DONE."

echo "Starting BAM files processing by prepare_fcnv_input.py"
for bam_file in $plasma_path/*.bam
do
    log_file=`echo $bam_file | sed -e 's/.bam/.log/g'`
    time -p $exec_path/prepare_fcnv_input.py $data_path/mp.phase.vcf $bam_file $data_path/__M.part.bam $data_path/__P.part.bam /dupa-filer/laci/centromeres $results_path > $log_file 2>&1 &
done
wait
echo "DONE."

#for bam_file in $plasma_path/*.bam
#do
#    log_file=`echo $bam_file | sed -e 's/.bam/.log/g'`
#    #qsub -q all.q -R y -pe parallel 6 -l h_vmem=10G -l h_rt=05:00:00 -S /bin/bash some_shell_script.sh
#    #qsub -q all.q -R y -pe parallel 6 -l h_vmem=10G -l h_rt=05:00:00 -S /usr/bin/python2 $exec_path/prepare_fcnv_input.py $data_path/mp.phase.vcf $bam_file $data_path/__M.part.bam $data_path/__P.part.bam /dupa-filer/laci/centromeres $results_path > $log_file 2>&1
#    
#done
#wait
