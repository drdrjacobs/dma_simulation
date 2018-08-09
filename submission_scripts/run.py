"""Submits run while varying parameter according to environment variable set 
by slurm jobarray."""

import numpy as np
import subprocess
import shlex
import os
import argparse
import pickle
import random
import time

# dimension 2 or 3
dimension = 2
# name of main parameter that is going to change from run to run
parameter = "fraction_max_kappa"
values = [0.1, 1.0]
# these params change along with the main parameter
params_dict = {}

executable_base = "ga_simulation"
id_key = "SLURM_ARRAY_TASK_ID"
# carry one decmial for 1.0e-3 parameter formatting
value_width = 1
# log format prefix$(id) in 010 format to log_width number of digits 
log_prefix = "log_"
log_width = 3
# folder is id_$(id)_$(parameter)_$(value), id 010 format above, 
# value is in 1.0e-6 format
params_out_file = "params.txt"

parser = argparse.ArgumentParser(description = __doc__)
parser.add_argument("executable_folder", 
                    help = "absolute path to executable folder")
parser.add_argument("params", help = "path to params.txt")
parser.add_argument("--global_rng_seed", type = int, help = "global rng " +
                    "seed added to id to get local rng seed, required " +
                    "for jobarray with multiple tasks")
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
value = values[id]
value_string = "{0:.{1}E}".format(value, value_width)
folder = "id_" + id_string + "_" + parameter + "_" + value_string
if not os.path.exists(folder):
    os.mkdir(folder)

# logfile
log_file = log_prefix + id_string
with open(folder + "/" + log_file, "w") as log:
    log.write("id: " + str(id) + "\n")
    log.write(parameter + ": " + str(value) + "\n")
    if args.global_rng_seed is not None:
        global_rng_seed = args.global_rng_seed
    elif id == 0:
        # either not run through job array or only has one task
        global_rng_seed = int(time.time())
        log.write("Setting global rng seed based on system time.")
    else:
        raise Exception("If global rng seed not set cannot run " +
                        "multiple job array tasks!")
    # all jobs have different seeds
    rng_seed = global_rng_seed + id
    log.write("global_rng_seed: " + str(global_rng_seed) + "\n")
    log.write("rng_seed: " + str(rng_seed) + "\n")

    # write modified params
    with open(folder + "/" + params_out_file, "w") as params_out:
        with open(args.params) as params_in:
            lines = [l.strip() for l in params_in]
            lines = [l for l in lines if l and l[0] != "#"]
            for l in lines:
                l = l.split()
                # format of line is type parameter = value
                parameter_type = l[0]
                p = l[1]
                v = l[3]
                if p == "seed":
                    v = rng_seed
                elif p == parameter:
                    v = value
                elif p in params_dict:
                    v = params_dict[p][id]
                new_l = "{} {} = {}\n".format(parameter_type, p, v)
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

