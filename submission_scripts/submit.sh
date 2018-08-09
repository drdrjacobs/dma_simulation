#!/bin/bash
export GLOBAL_RNG_SEED=$(date +%s)
values="$(grep 'values =' run.py)"
# count commas to get number of values - 1
count="$(echo $values | grep -o , | wc -l)"
# submit number of values jobs
sbatch --array 0-$count run.sh
