#!/bin/bash

time ./sim_duplicate.sh 10000000 100000 P &
time ./sim_duplicate.sh 10000000 1000000 P &
time ./sim_duplicate.sh 10000000 10000000 P &

time ./sim_duplicate.sh 10000000 100000 M &
time ./sim_duplicate.sh 10000000 1000000 M &
time ./sim_duplicate.sh 10000000 10000000 M &

