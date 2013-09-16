#!/bin/bash
for i in 1
do
  pypy ../data_generator/data_generator2.py $i &
done
#pypy ../data_generator/data_generator2.py $1 &
