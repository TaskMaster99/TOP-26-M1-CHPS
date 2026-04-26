#!/bin/bash

sudo sh -c 'echo -1 > /proc/sys/kernel/perf_event_paranoid'

mode=("STAT" "RECORD" "REPORT" "HOTSPOT")
size=$2


if [[ ${mode[0]} == ${1} ]]
then
    perf stat -e branch-instructions,\
    cpu-cycles,\
    dTLB-loads,\
    dTLB-load-misses,\
    iTLB-loads,\
    iTLB-load-misses,\
    L1-icache-load-misses,\
    L1-icache-loads,\
    branch-misses,\
    branch-load-misses,\
    cpu-clock,\
    L1-dcache-loads,\
    L1-dcache-load-misses,\
    L1-dcache-prefetches \
    -- \
    mpirun -np ${size} ./build/top.lbm-exe config.txt
elif [[ ${mode[1]} == ${1} ]]
then
    rm -r profiling/data/perf*
    mkdir -p $(pwd)/profiling/data
    mpirun -np ${size} "$(pwd)/profiling/script/sample.sh"
elif [[ ${mode[2]} == ${1} ]]
then
    perf report -i "$(pwd)/profiling/data/perf.${2}.data"
else
    sudo ./profiling/bin/hotspot.AppImage "$(pwd)/profiling/data/perf.${2}.data"
fi
