#!/bin/bash
#PBS -N bve_time_9
#PBS -A UMIC0093
#PBS -l walltime=12:00:00
#PBS -q regular
#PBS -j oe
#PBS -k eod
#PBS -m abe
#PBS -M cygnari@umich.edu
#PBS -l select=1:ncpus=36:mpiprocs=36
#PBS -l place=group=rack

export TMPDIR=/glade/scratch/$USER/temp
mkdir -p $TMPDIR

mpirun -np 2 ./driver > $TMPDIR/run_out9.txt
