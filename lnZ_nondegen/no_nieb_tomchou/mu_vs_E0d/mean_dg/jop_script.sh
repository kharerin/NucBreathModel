#!/bin/bash
##PBS -j oe
##PBS -q serial
##PBS -N ms0p1x
##PBS -l walltime=10:00:00

##echo "Running on: " 
##cat ${PBS_NODEFILE}

##echo "Program Output begins: "

##cd ${PBS_O_WORKDIR}

gcc hung_soft_lnZ_set.c -lm -o den.out
nohup ./den.out > denout.txt &
