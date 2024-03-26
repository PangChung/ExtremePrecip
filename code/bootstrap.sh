#!/bin/bash

#PBS -l ncpus=16
#PBS -l mem=124gb
#PBS -l walltime=12:00:00
#PBS -j oe 
##PBS -J 1-2400
#PBS -m a
#PBS -M peng.zhong@unsw.edu.au 

# cd ${PBS_O_WORKDIR}

cd ~/R/ExtremePrecip/
module load r/4.3.1

# region=$(( PBS_ARRAY_INDEX / 301 + 1 ))
# boot_ind=$(( PBS_ARRAY_INDEX % 301 + 1))

# region=$(( PBS_ARRAY_INDEX / 8 + 1 ))
# boot_ind=$(( PBS_ARRAY_INDEX % 8 + 1))

# Rscript code/bootstrap.R "idx.region=${region};bootstrap.ind=${boot_ind};computer=\"hpc\"" 

# # Array to store the process IDs of each job
declare -a pids

# date
# Maximum number of parallel jobs
max_jobs=4

# Function to run a job with nohup
run_job() {
    local region=$1
    local boot_ind=$2
    local output_file="output_${region}_${boot_ind}.txt"
    nohup Rscript code/bootstrap.R "idx.region=${region};bootstrap.ind=${boot_ind};computer=\"hpc\""  > "$output_file" 2>&1  &
    echo $!
    pids+=($!)
}

is_running() {
    local pid=$1
    kill -0 "$pid" >/dev/null 2>&1
}

# Run multiple jobs in parallel
for region in {1..8}; do
    while [ ${#pids[@]} -ge $max_jobs ]; do
        sleep 1
        # Remove finished jobs from the pids array
        for ((i=0; i<${#pids[@]}; i++)); do
            if ! is_running "${pids[$i]}"; then
                unset 'pids[$i]'
            fi
        done
        pids=("${pids[@]}")    
    done
    echo "Running job $region"
    run_job $region $boot
done

date
