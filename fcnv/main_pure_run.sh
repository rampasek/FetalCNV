#!/bin/bash

pipeline_stage=$1
chr=$2
ff=$3

#plasma_path="/dupa-filer/laci/I1/pure/"
#results_path="/dupa-filer/laci/I1/pure/ff$ff/"
results_path="/dupa-filer/laci/I1/pureAll/ff$ff/"
exec_path="/dupa-filer/laci/bin/"


job_count=0
flag=0
for x in 22 X {1..21} #3 4 8  #{1..2} {5..7} #{9..22}
do
#for bam_file in $plasma_path/*.bam
#do
    chr=chr$x
    data_path="/dupa-filer/laci/I1/chr$x/"
    #bam_file=$data_path/__plasma-chr$x-0.$ff'rate.sort.bam'
    bam_file=$data_path/__plasma.part.bam
    log_file=$results_path/log_ff$ff.chr$x.log
    
    if [[ "$flag" == 10 ]] ; then
        flag=0
        sleep 500
    fi
    flag=$(($flag+1))
    
    #echo $log_file
    echo "qsub -q all.q -V -R y -pe parallel 8 -l h_vmem=6G -l h_rt=12:00:00 -o $log_file -e $log_file -S /usr/bin/pypy $exec_path/prepare_pure_input.py $data_path/mp.phase.vcf $bam_file $data_path/__M.part.bam $data_path/__P.part.bam $data_path/regions.bed /dupa-filer/laci/centromeres $results_path"
    #qsub -q all.q -V -R y -pe parallel 8 -l h_vmem=6G -l h_rt=12:00:00 -o $log_file -e $log_file -S /usr/bin/pypy $exec_path/prepare_pure_input.py $data_path/mp.phase.vcf $bam_file $data_path/__M.part.bam $data_path/__P.part.bam $data_path/regions.bed /dupa-filer/laci/centromeres $results_path
    
    qsub -q all.q -V -R y -pe parallel 8 -l h_vmem=6G -l h_rt=12:00:00 -o $log_file -e $log_file -S /usr/bin/pypy $exec_path/prepare_pureAll_input.py $data_path/mp.phase.vcf $bam_file $data_path/__M.part.bam $data_path/__P.part.bam $data_path/regions.bed /dupa-filer/laci/centromeres $results_path
    
    sleep 10
    
done
echo "DONE."
