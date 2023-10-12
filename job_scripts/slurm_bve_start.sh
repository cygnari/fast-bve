#!/bin/bash
# (See https://arc-ts.umich.edu/greatlakes/user-guide/ for command details)

# Set up batch job settings
#SBATCH --job-name=bve
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=36
#SBATCH --exclusive
#SBATCH --time=12:00:00
#SBATCH --account=krasny0
#SBATCH --partition=standard

export TMPDIR=/home/$USER/temp
mkdir -p $TMPDIR

mpirun -np 36 --bind-to core:overload-allowed ../build/executables/driver > $TMPDIR/run_out.txt
