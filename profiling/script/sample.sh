#!/bin/bash

rank=${PMI_RANK}
perf record --call-graph dwarf -o $(pwd)/profiling/data/perf.${rank}.data -g ./build/top.lbm-exe config.txt
