 #!/bin/bash
export OMP_PROC_BIND=spread
export OMP_PLACES=threads
export OMP_NUM_THREADS=8
size=$1

mpirun -np ${size} --bind-to core ./build/top.lbm-exe config.txt