#!/bin/bash
# check the number of free cores and run program when enough cores are free

# Desired number of free cores
REQUIRED_FREE_CORES=16

# Total cores on this computer
TOTAL_CORES=32

# The program to run
PROGRAM="./rebound &"

# Interval between checks in seconds
CHECK_INTERVAL=10

# Function to get the number of free cores
get_free_cores() {

  # Run top in batch mode and capture the running cores
  top_output=$(top -b | head -n 5)
  RUNNING_CORES=$(echo "$top_output" | grep -oP '\d+(?= running)')

  # Calculate the number of free cores
  FREE_CORES=$((TOTAL_CORES - RUNNING_CORES))

  echo $FREE_CORES
}

# Main loop
while true; do
  # Get the number of free cores
  FRE_CORES=$(get_free_cores)

  # Check if the number of free cores meets the requirement
  if [[ FRE_CORES -ge REQUIRED_FREE_CORES ]]; then
    echo "Starting the program as $FRE_CORES cores are free."
    $PROGRAM &
    exit 0
  else
    echo "$FRE_CORES cores are free. Waiting for at least $REQUIRED_FREE_CORES free cores."
  fi

  # Wait before the next check
  sleep $CHECK_INTERVAL
done

