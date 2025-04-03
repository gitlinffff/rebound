#!/bin/bash

# List of r_dust values to loop over
r_dust=(1.00000e-04 1.41254e-04 1.99526e-04 2.81838e-04 3.98107e-04 \
	5.62341e-04 7.94328e-04 1.12202e-03 1.58489e-03 2.23872e-03 \
	3.16228e-03 4.46684e-03 6.30957e-03 8.91251e-03 1.25893e-02 \
	1.77828e-02 2.51189e-02 3.54813e-02 5.01187e-02 7.07946e-02 \
	1.00000e-01)

ROOTDIR=$(pwd)
SUBTEMP="submit.bro-28"
EXE="rebound"
LIB="librebound.so"

i=1
# Loop and create PBS scripts
for r in "${r_dust[@]}"; do
    tag=$(echo $r | sed 's/e/E/')  # Replace 'e' with 'E' for job name clarity
    script="submit_r${tag}.pbs"
    run_id=$(printf "%03d" $i)
    rundir="${ROOTDIR}/run_${run_id}"
    mkdir -p "${rundir}"

    # Replace placeholder with actual value
    cd "$rundir"
    sed -e "s/{{RDUST}}/${r}/g" \
	-e "s/{{JOBID}}/radius${run_id}/g" \
	"${ROOTDIR}/${SUBTEMP}" > "$script"
    
    ln -s "${ROOTDIR}/${EXE}"
    ln -s "${ROOTDIR}/${LIB}"

    # Submit the job
    qsub "$script"
    
    cd "$ROOTDIR"
    ((i++))
done

