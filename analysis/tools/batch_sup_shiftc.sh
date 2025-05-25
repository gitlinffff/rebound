#!/bin/bash

ROOTDIR=$(pwd)

# Run this on your **local machine**
for i in $(seq 14 16); do   # directory run_002 to run_021
    
    id=$(printf "%03d" $i)
    
    local_dir="${ROOTDIR}/run_${id}"
    remote_file="lfe:/u/lli22/ejecta_size_exp/run_${id}/particles.txt"
    #remote_file="lfe:/u/lli22/test/run_${id}/particles.txt"

    echo ${local_dir}

    mkdir -p "${local_dir}"
    cd "${local_dir}"

    sup shiftc ${remote_file} .
    
    cd "$ROOTDIR"
    
    sleep 5
done
