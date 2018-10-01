#!/bin/bash
export GLOBAL_RNG_SEED=$(date +%s)
# python script writes file count.sh
# this file contains the variable COUNT that is the number of jobs to run
python set_up_job_input.py
source count.sh
# submit number of values jobs
sbatch --array 0-$COUNT run.sh
/bin/rm count.sh
