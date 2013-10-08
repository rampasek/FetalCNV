#!/bin/bash

pipeline_stage=$1
chr=$2
data_path="/dupa-filer/laci/I1/$chr/"
plasma_path="/dupa-filer/laci/I1/$chr/sim_plasma/"
results_path="/dupa-filer/laci/I1/$chr/fcnv_data/"
exec_path="/dupa-filer/laci/bin/"

if [[ ! -d "$data_path" || ! -d "$plasma_path" ]]; then
    echo "ERROR: data_path or plasma_path doesn't exist!"
    exit
fi

length_set='10000000 1000000 500000 100000'
haplo_set='A B'
src_set='P M'
#begin_set='10000000 35000000 47000000'
begin_set='10000000 75000000 150000000 190000000'

export SGE_O_PATH=$PATH

if [[ "$pipeline_stage" == 1 ]] ; then
echo "Starting CNV plasma files simulations"
for begin in $begin_set
do
	for length in $length_set
	do
		for haplo in $haplo_set
		do
			qsub -q all.q -V -R y -l h_vmem=5G -l h_rt=05:00:00 -S /bin/bash $exec_path/sim_deletion.sh $chr $begin $length $haplo $data_path $plasma_path $exec_path
			echo $begin' '$length' '$haplo

			for src in $src_set
			do
				qsub -q all.q -V -R y -l h_vmem=5G -l h_rt=05:00:00 -S /bin/bash $exec_path/sim_duplicate.sh $chr $begin $length $src $haplo $data_path $plasma_path $exec_path
				echo $begin' '$length' '$src' '$haplo
			done
		done
	done
done
echo "DONE."
exit
fi

#$exec_path/sim_duplicate.sh 10000000 100000 P B $data_path $plasma_path $exec_path &
#$exec_path/sim_deletion.sh 10000000 100000 A $data_path $plasma_path $exec_path &
#wait

if [[ ! -d "$results_path" ]]; then
    echo "ERROR: results_path doesn't exist!"
    exit
fi

if [[ "$pipeline_stage" == 2 ]] ; then
echo "Starting BAM files processing by prepare_fcnv_input.py"
for bam_file in $plasma_path/*.bam
do
    log_file=`echo $bam_file | sed -e 's/.bam/.log/g' | sed -e 's/:/-/g'`
    qsub -q all.q -V -R y -pe parallel 7 -l h_vmem=5G -l h_rt=05:00:00 -o $log_file -e $log_file -S /usr/bin/python2 $exec_path/prepare_fcnv_input.py $data_path/mp.phase.vcf $bam_file $data_path/__M.part.bam $data_path/__P.part.bam $data_path/regions.bed /dupa-filer/laci/centromeres $results_path
    #time -p $exec_path/prepare_fcnv_input.py $data_path/mp.phase.vcf $bam_file $data_path/__M.part.bam $data_path/__P.part.bam /dupa-filer/laci/centromeres $results_path > $log_file 2>&1 &
done
wait
echo "DONE."
exit
fi


if [[ "$pipeline_stage" == 3 ]] ; then
echo "Running delOrigin.py"
for del_tgt_file in $results_path/*delete.target.txt
do
    delOrigin.py $data_path/trio.phase.vcf $del_tgt_file &
done
wait
echo "DONE."
fi

