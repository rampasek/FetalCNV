#!/bin/bash

pypy data_generator.py T$1 &
pypy data_generator.py $1 &
