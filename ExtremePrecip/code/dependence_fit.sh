#!/bin/bash
cd ~/Desktop/ExtremePrecip/

# Array to store the process IDs of each job
declare -a pids

# Maximum number of parallel jobs
max_jobs=5

# Function to run a job with nohup
run_job() {
    local region=$1
    local season=$2
    local norm=$3
    local output_file="output_${region}_${season}_${norm}.txt"
    #local error_file="error_${region}_${season}_${norm}.txt"
    nohup Rscript code/dependence_fit.R "idx.region=$region;season.ind=$season;norm.ind=$norm" > "$output_file" 2>&1  &
    echo $!
    pids+=($!)
}

is_running() {
    local pid=$1
    kill -0 "$pid" >/dev/null 2>&1
}

# Run multiple jobs in parallel
for region in {1..8}; do
    for season in {1..4}; do
        for norm in {1..2}; do
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
            echo "Running job $region $season $norm"
            run_job $region $season $norm
        done
    done
done