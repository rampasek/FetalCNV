#!/bin/bash

data_path='/dupa-filer/laci/I1/chr20/'
plasma_path='/dupa-filer/laci/I1/chr20/sim_plasma/'
results_path='/dupa-filer/laci/I1/chr20/fcnv_data/'
exec_path='/dupa-filer/laci/bin/'

#time ./sim_duplicate.sh 10000000 100000 P &
$exec_path/sim_duplicate.sh 10000000 100000 P B $data_path $plasma_path $exec_path 
#$exec_path/sim_deletion.sh 10000000 1000000 A $data_path $plasma_path $exec_path

bam_file=$plasma_path/*.bam
log_file=`echo $bam_file | sed -e 's/.bam/.log/g'`
time -p $exec_path/prepare_fcnv_input.py $data_path/mp.phase.vcf $bam_file $data_path/__M.part.bam $data_path/__P.part.bam /dupa-filer/laci/centromeres $results_path > $log_file 2>&1 &

#for bam_file in $plasma_path/*.bam
#do
#    log_file=`echo $bam_file | sed -e 's/.bam/.log/g'`
#    #qsub -q all.q -R y -pe parallel 6 -l h_vmem=10G -l h_rt=05:00:00 -S /bin/bash some_shell_script.sh
#    #qsub -q all.q -R y -pe parallel 6 -l h_vmem=10G -l h_rt=05:00:00 -S /usr/bin/python2 $exec_path/prepare_fcnv_input.py $data_path/mp.phase.vcf $bam_file $data_path/__M.part.bam $data_path/__P.part.bam /dupa-filer/laci/centromeres $results_path > $log_file 2>&1
#    
#done
#wait
