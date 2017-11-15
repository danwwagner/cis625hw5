#!/bin/bash
#$ -l killable
#$ -cwd
#$ -j y
#$ -l h_rt=01:00:00
#$ -l mem=1G
#$ -N Q2_MPI_8
#$ -m as
#$ -M danwagner@ksu.edu

mpirun ./q2-hybrid $1 $2 $3
