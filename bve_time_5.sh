#!/bin/bash
#PBS -N bve_time_5
#PBS -A UMIC0093
#PBS -l walltime=12:00:00
#PBS -q main
#PBS -j oe
#PBS -k eod
#PBS -m abe
#PBS -M cygnari@umich.edu
#PBS -l select=8:ncpus=64:mpiprocs=64
#PBS -l place=group=rack

export TMPDIR=/glade/derecho/scratch/$USER/temp
mkdir -p $TMPDIR

mpirun -np 512 ./driver > $TMPDIR/run_out5.txt
