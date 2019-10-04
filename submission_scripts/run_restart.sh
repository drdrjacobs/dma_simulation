#!/bin/bash
#SBATCH -J dma # name 
#SBATCH -p CLUSTER,debug # partition
#SBATCH -t 01-00:00 # time (D-HH:MM)
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=63000 # in mb, must have lots of memory for large clusters
# must have .git access
executable_folder="/home/drjacobs/dendrites/dma_simulation/bin" 
current=$(pwd)
cd $executable_folder
git log -1 HEAD --decorate > $current/version.txt
cd $current
python3 run_restart.py $executable_folder
