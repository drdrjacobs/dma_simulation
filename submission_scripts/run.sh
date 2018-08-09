#!/bin/bash
#SBATCH -J ga_simulation # name 
#SBATCH -p any # partition
#SBATCH -t 0-00:30 # time (D-HH:MM)
#SBATCH --mem=10000 # in mb, must have lots of memory for large clusters
# must have .git access
executable_folder="/home/djacobso/dendrites/ga_simulation/bin" 
current=$(pwd)
cd $executable_folder
git log -1 HEAD --decorate > $current/version.txt
cd $current
python3 run.py $executable_folder params.txt \
    --global_rng_seed $GLOBAL_RNG_SEED

