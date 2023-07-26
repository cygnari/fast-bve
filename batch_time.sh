#!/bin/bash
export MPI_RANKS=1
export NODES=1
qsub bve_time_0.sh
qsub bve_time_1.sh
qsub bve_time_2.sh
qsub bve_time_3.sh
qsub bve_time_4.sh
qsub bve_time_5.sh
qsub bve_time_6.sh
qsub bve_time_7.sh
qsub bve_time_8.sh
qsub bve_time_9.sh
