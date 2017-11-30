#!/bin/bash
#----------------------------------------------------
#SBATCH -J myjob           # Job name
#SBATCH -o myjob.o%j       # Name of stdout output file
#SBATCH -e myjob.e%j       # Name of stderr error file
#SBATCH -p normal          # Queue (partition) name
#SBATCH -N 4               # Total # of nodes
#SBATCH -n 32              # Total # of mpi tasks
#SBATCH -t 01:30:00        # Run time (hh:mm:ss)

# Launch MPI code...


for i in 1 2 4 8 16 32 64 128 256;do for t in 1 2 4 8 16 32 64 128 256;do n=$( printf "%.0f" $(echo "sqrt ( $t ) * 1000" | bc -l) );ibrun -n $i ./mpi_nbody3 -n $n 2>/dev/null;done > $i-weak-output.log;done;
