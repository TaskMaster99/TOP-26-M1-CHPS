#!/bin/bash



# Configure
cmake -B build

# Compile
cmake --build build -t top.lbm-exe
cmake --build build -t top.display
