#!/bin/bash

#PBS -l select=16:ncpus=1:mem=124gb
#PBS -l walltime=11:59:59 ## ideal
#PBS -j oe ##merge output and error into a single file
#PBS -J 1-300 ## array jobs 
#PBS -M peng.zhong@unsw.edu.au ## notify me via email when job ends (a) or terminated because of error (e)
#PBS -m ae ## 

# cd ${PBS_O_WORKDIR}

# ./myprogram ${PBS_ARRAY_INDEX}.dat

cd ~/Desktop/ExtremePrecip/

# Array to store the process IDs of each job
declare -a pids

date
# Maximum number of parallel jobs
max_jobs=4

# Function to run a job with nohup
run_job() {
    local region=$1
    local boot_ind=$2
    local output_file="output_${region}_${boot_ind}.txt"
    nohup Rscript code/bootstrap.R "idx.region=$region;bootstrap.ind=$boot_ind" > "$output_file" 2>&1  &
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
    run_job $region $1
done

date
