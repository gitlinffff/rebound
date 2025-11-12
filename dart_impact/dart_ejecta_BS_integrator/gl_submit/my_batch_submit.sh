#!/bin/bash

# List of r_dust values to loop over
# 1e-6 ~ 1.84206997e-03  33 radius
r_dust=(1.00000000e-06 1.26485522e-06 1.59985872e-06 2.02358965e-06 \
	2.55954792e-06 3.23745754e-06 4.09491506e-06 5.17947468e-06 \
	6.55128557e-06 8.28642773e-06 1.04811313e-05 1.32571137e-05 \
	1.67683294e-05 2.12095089e-05 2.68269580e-05 3.39322177e-05 \
	4.29193426e-05 5.42867544e-05 6.86648845e-05 8.68511374e-05 \
	1.09854114e-04 1.38949549e-04 1.75751062e-04 2.22299648e-04 \
	2.81176870e-04 3.55648031e-04 4.49843267e-04)

ROOTDIR=$(pwd)
SUBTEMP="submit-reb.sub"
EXE="rebound"
LIB="librebound.so"

i=1
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
