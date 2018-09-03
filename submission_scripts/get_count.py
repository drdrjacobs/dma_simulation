"""Finds number of jobs to run based on run.py and writes this to a file that
bash can read in."""

with open("run.py") as f:
    lines = [_ for _ in f]
values_line = [_ for _ in lines if _[:9] == "values = "]
values_line = values_line[0]
repeats_line = [_ for _ in lines if _[:9] == "repeats = "]
repeats_line = repeats_line[0]
# syntax is "values = list(np.logspace(min, max, N))"
count = int(values_line.split(",")[-1].split(")")[0])
# syntaxs is "repeats = N" 
count = count * int(repeats_line.split("=")[-1])
# jobs start at index 0
count -= 1
# write out to bash file
with open("count.sh", "w") as f:
    f.write("COUNT={}\n".format(count))

    
