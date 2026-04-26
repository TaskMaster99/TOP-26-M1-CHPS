 #!/bin/bash

size=$1

mpirun -np ${size} ./build/top.lbm-exe config.txt