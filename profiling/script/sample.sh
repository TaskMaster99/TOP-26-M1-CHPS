#!/bin/bash

rank=${PMI_RANK}
perf record --call-graph dwarf -o $(pwd)/profiling/data/perf.${rank}.data -g timeout 1m ./build/top.lbm-exe config.txt
