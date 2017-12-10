#!/bin/bash
## sbatch -J unboxzip -p normal -N 1 -n 32 -t 05:00:00 unbox.sh

for i in $(ls $WORK/data/raw/); do bunzip2 $WORK/data/raw/$i;done
