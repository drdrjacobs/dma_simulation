#!/bin/bash
# python script writes file count.sh
# this file contains the variable COUNT that is the number of jobs to run
python set_up_restart.py
source count.sh
# submit number of values jobs
sbatch --array 0-$COUNT%10 run_restart.sh
/bin/rm count.sh
