#!/bin/bash
#PBS -N bve_run
#PBS -A UMIC0093
#PBS -l walltime=12:00:00
#PBS -q main
#PBS -j oe
#PBS -k eod
#PBS -m abe
#PBS -M cygnari@umich.edu
#PBS -l select=1:ncpus=1:mpiprocs=1
#PBS -l place=group=rack

module load linaro-forge/23.0

export TMPDIR=/glade/derecho/scratch/$USER/temp
mkdir -p $TMPDIR

perf-report --mpi -n 1 ../build/executables/single_rhs
