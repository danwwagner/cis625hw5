#!/bin/bash
# Varying the number of processes, up to 16

./mass_qsub.sh 1 $1 $2
./mass_qsub.sh 2 $1 $2 
./mass_qsub.sh 4 $1 $2 
./mass_qsub.sh 8 $1 $2 
./mass_qsub.sh 16 $1 $2

