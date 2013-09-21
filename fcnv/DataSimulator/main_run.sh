$data_path
$results_path
$exec_path

time ./sim_duplicate.sh 10000000 100000 P &
$exec_path/duplication/sim_duplicate.sh 10000000 100000 P $data_path $results_path &
#$exec_path/deletion/sim_deletion.sh 10000000 1000000 $data_path $results_path &

