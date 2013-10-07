#!/bin/bash
data_path='/dupa-filer/laci/I1/chr20/'
plasma_path='/dupa-filer/laci/I1/chr20/test/'
results_path='/dupa-filer/laci/I1/chr20/test/'
exec_path='/dupa-filer/laci/bin/'

time $exec_path/sim_deletion.sh 35000000 100000 P $data_path $plasma_path $exec_path

#time ./sim_duplicate.sh 10000000 100000 P &
#time ./sim_duplicate.sh 10000000 1000000 P &
#time ./sim_duplicate.sh 10000000 10000000 P &

#time ./sim_duplicate.sh 10000000 100000 M &
#time ./sim_duplicate.sh 10000000 1000000 M &
#time ./sim_duplicate.sh 10000000 10000000 M &

