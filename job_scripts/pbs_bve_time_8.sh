#!/bin/bash
#PBS -N bve_time_8
#PBS -A UMIC0093
#PBS -l walltime=12:00:00
#PBS -q main
#PBS -j oe
#PBS -k eod
#PBS -m abe
#PBS -M cygnari@umich.edu
#PBS -l select=1:ncpus=8:mpiprocs=8
#PBS -l place=group=rack

export TMPDIR=/glade/derecho/scratch/$USER/temp
mkdir -p $TMPDIR

mpiexec ../build/executables/single_rhs > $TMPDIR/run_out108.txt
