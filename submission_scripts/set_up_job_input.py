"""Sets up inputs for jobs. Specifies what parameters each job will have."""

import numpy as np
import pickle as pkl

# dimension 2 or 3 
dimension = 2
# name of main parameter that is going to change from run to run
parameter = "kappa"
values = np.array([25, 250, 2500, 16 * 2500])
values = np.sqrt(values)
# number of repeats of each value to run
repeats = 100
# these params change along with the main parameter
params_dict = {}

# sort everything appropriately
values = list(values) * repeats
for k, v in params_dict.items():
    v = v * repeats
    params_dict[k] = [_ for value, _ in sorted(zip(values, v))]
values = sorted(values)

out_dict = {"parameter": parameter, "values": values, 
            "params_dict": params_dict, "dimension": dimension}
# pickle parameters for later use within each job
with open("values.pkl", "wb") as f:
    pkl.dump(out_dict, f)

# write this out to specify how many jobs are being run in total
count = len(values) - 1
with open("count.sh", "w") as f:
    f.write("COUNT={}\n".format(count))

    
