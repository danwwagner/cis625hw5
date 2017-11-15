#!/bin/bash
#$ -l killable
#$ -j y
#$ -l mem=1G
#$ -m abe
#$ -M danwagner@ksu.edu

for i in {1}
	do
		qsub -pe mpi-fill $1 -q \*@@elves mpi_qsub.sh $2 $3 8
	done
