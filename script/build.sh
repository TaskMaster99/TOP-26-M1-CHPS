#!/bin/bash

export OMP_PROC_BIND=spread
export OMP_PLACES=threads


# Configure
cmake -B build

# Compile
cmake --build build -t top.lbm-exe
cmake --build build -t top.display
