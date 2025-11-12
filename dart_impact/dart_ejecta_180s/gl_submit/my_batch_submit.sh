#!/bin/bash

# List of r_dust values to loop over
# 1.20679264e-02 ~ 1e-01  10 radius
r_dust=(1.20679264e-02 1.52641797e-02 1.93069773e-02 2.44205309e-02 \
	3.08884360e-02 3.90693994e-02 4.94171336e-02 6.25055193e-02 \
	7.90604321e-02 1.00000000e-01)

ROOTDIR=$(pwd)
SUBTEMP="submit-reb.sub"
EXE="rebound"
LIB="librebound.so"

i=41  # run index from which to start
# Loop and create PBS scripts
for r in "${r_dust[@]}"; do
	tag=$(echo $r | sed 's/e/E/')  # Replace 'e' with 'E' for job name clarity
	run_id=$(printf "%03d" $i)
	script="submit_${run_id}.sub"
	rundir="${ROOTDIR}/run_${run_id}"
	mkdir -p "${rundir}"

	# Replace placeholder with actual value
	cd "$rundir"
	sed -e "s/{{RDUST}}/${r}/g" \
	-e "s/{{JOBID}}/radius_${run_id}/g" \
	"${ROOTDIR}/${SUBTEMP}" > "$script"

	ln -s "${ROOTDIR}/${EXE}"
	ln -s "${ROOTDIR}/${LIB}"

	# Submit the job
	sbatch "$script"

	cd "$ROOTDIR"
	((i++))
done
