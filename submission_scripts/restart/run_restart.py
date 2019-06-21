"""Submits restarts with potentially changed parameters according to 
environment variable set by slurm jobarray."""

import numpy as np
import subprocess
import shlex
import os
import argparse
import pickle
import random
import time
import IPython as ipy
import pickle as pkl
import glob

with open("restart_values.pkl", "rb") as f:
    data = pkl.load(f)
# dimension 2 or 3
dimension = data["dimension"]
# these params change in the restarted run
params_dict = data["params_dict"]

executable_base = "ga_simulation"
id_key = "SLURM_ARRAY_TASK_ID"
# log format prefix$(id) in 010 format to log_width number of digits 
log_prefix = "restart_log_"
log_width = 3
params_file = "params.txt"

parser = argparse.ArgumentParser(description = __doc__)
parser.add_argument("executable_folder", 
                    help = "absolute path to executable folder")
args = parser.parse_args()

if args.executable_folder[-1] != "/":
    args.executable_folder += "/"

# get jobarray id
if id_key in os.environ:
    id = int(os.environ[id_key])
else:
    # not being run by job array
    id = 0
id_string = "{0:0{1}}".format(id, log_width)

folders = sorted(glob.glob("id_*"))
folder = folders[id]

# logfile
log_file = log_prefix + id_string
with open(folder + "/" + log_file, "w") as log:
    log.write("id: " + str(id) + "\n")
    # sort potential restart files
    restart_files = glob.glob(folder + "/restart_*.ser")
    restart_files = [_.split("/")[-1] for _ in restart_files]
    restart_orderings = [int(_.split("_")[-1][:-4]) for _ in restart_files] 
    restart_files = [restart_file for _, restart_file in 
                     sorted(zip(restart_orderings, restart_files))]
    # need to pick restart file, start from second to newest file in case 
    # newest did not finish writing out
    if len(restart_files) >= 2:
        restart_file = restart_files[-2]
    else:
        restart_file = restart_files[0]
    # write modified params
    os.rename(folder + "/" + params_file, folder + "/old_" + params_file)
    with open(folder + "/" + params_file, "w") as params_out:
        with open(folder + "/old_" + params_file) as params_in:
            lines = [l.strip() for l in params_in]
            lines = [l for l in lines if l and l[0] != "#"]
            has_restart = False
            for l in lines:
                l = l.split()
                # format of line is type parameter = value
                parameter_type = l[0]
                p = l[1]
                v = l[3]
                if p == "restart_path":
                    has_restart = True
                    v = restart_file
                elif p in params_dict:
                    v = params_dict[p][id]
                new_l = "{} {} = {}\n".format(parameter_type, p, v)
                params_out.write(new_l)
            if not has_restart:
                new_l = "std::string restart_path = {}\n".format(restart_file)
                params_out.write(new_l)

    # os call to run
    os.chdir(folder)
    call = ("time " + args.executable_folder + executable_base + 
            "_{}d".format(dimension))
    log.write("call: " + call + "\n")
    log.write("Done with Python.\n\n")
    call = shlex.split(call)
    
# and done!
log = open(log_file, "a")
subprocess.call(call, stdout = log)

