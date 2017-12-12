#!/bin/bash
#----------------------------------------------------
#SBATCH -J myjob           # Job name
#SBATCH -o myjob.o%j       # Name of stdout output file
#SBATCH -e myjob.e%j       # Name of stderr error file
#SBATCH -p normal          # Queue (partition) name
#SBATCH -N 4               # Total # of nodes
#SBATCH -n 32              # Total # of mpi tasks
#SBATCH -t 01:30:00        # Run time (hh:mm:ss)

# Launch OpenMP code benchamark...
# sbatch -J volume -p normal -N 1 -n 68 -t 05:00:00 benchmark.sh


for t in 1 2 4 6 8 16 32 64; do export OMP_NUM_THREADS=$t; for i in $(ls $WORK/data/05/);do $HOME/fixanalyzer -p $WORK/data/05/$i -t 52 -n 4 -m 8 -b > $HOME/output/$t-${i:11:17}-buffer.log; done; done
