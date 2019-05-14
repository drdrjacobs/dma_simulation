"""Specifies what parameters will change for restarted jobs if any."""

import glob
import pickle as pkl

# dimension 2 or 3 
dimension = 2
count = len(glob.glob("id_*"))
params_dict = {"cluster_size": [300e6] * count}

out_dict = {"params_dict": params_dict,
            "dimension": dimension}

# pickle parameters for later use within each job
with open("restart_values.pkl", "wb") as f:
    pkl.dump(out_dict, f)
# write out how many jobs will run
with open("count.sh", "w") as f:
    # subtract 1 since jobarray starts at zero
    f.write("COUNT={}\n".format(count - 1))
