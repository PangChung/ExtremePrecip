#!/bin/bash
cd ~/Desktop/ExtremePrecip/

# Array to store the process IDs of each job
declare -a pids

# Maximum number of parallel jobs
max_jobs=8

# Function to run a job with nohup
run_job() {
    local region=$1
    local season=$2
    local norm=$3
    nohup Rscript code/dependence_fit.R "idx.region=$region;season.ind=$season;norm.ind=$3" &
    pids+=($!)
}

# Run multiple jobs in parallel
for region in {1..8}; do
    for season in {1..4}; do
        for norm in {1..2}; do
            # Check if the maximum number of parallel jobs has been reached
            while [ ${#pids[@]} -ge $max_jobs ]; do
                sleep 1
            done
            run_job $region $season $norm
        done
    done
done

# Wait for all jobs to finish
for pid in "${pids[@]}"; do
    wait "$pid"
done

# Check if all jobs have finished
if [ ${#pids[@]} -eq N ]; then
    echo "All jobs have finished."
else
    echo "Some jobs are still running."
fi

echo "Your program has been finished" | mail -s "Your program has been finished" -a "From:zhongprw@gmail.com" -s smtp="smtp.gmail.com:587" -s smtp-use-starttls -s smtp-auth=login -s smtp-auth-user="zhongprw@gmail.com" -s "Marginal Fit" peng.zhong@unsw.edu.au 