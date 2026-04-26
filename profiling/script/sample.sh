#!/bin/bash

rank=${PMI_RANK}
perf record -o $(pwd)/profiling/data/perf.${rank}.data -g timeout 1m ./build/top.lbm-exe config.txt
