#!/bin/bash

/usr/bin/time -p pypy fcnv.py $2 > $1 2>> $1
