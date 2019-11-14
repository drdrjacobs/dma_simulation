'''This script is used to verify that various distributions are correct. This
is mostly done visually through plots.'''

import matplotlib.pyplot as plt
import numpy as np
import IPython as ipy
import glob
import scipy.stats

file = "generate_point_on_sphere_2d.txt"
data = np.loadtxt(file)
# y coord is first
angles = [np.arctan2(_[1], _[0]) for _ in data] 
plt.hist(angles, density = True, bins = 50, range = (-np.pi, np.pi), 
         label = "$N = " + str(len(angles)) + "$")
plt.plot([-np.pi, np.pi], np.array([1, 1]) / (2 * np.pi), label = "Uniform")
plt.title(file.replace("_", "\_").replace(".txt", ""))
plt.xlabel(r"$\theta$")
plt.ylabel(r"$P(\theta)$")
plt.legend()
plt.show()

file = "generate_point_on_sphere_3d.txt"
data = np.loadtxt(file)
epsilon = 1e-8
# y coord is first
phis = [np.arctan2(_[1], _[0]) for _ in data if _[0]**2 + _[1]**2 > epsilon] 
plt.hist(phis, density = True, bins = 50, range = (-np.pi, np.pi),
         label = "$N = " + str(len(phis)) + "$")
plt.plot([-np.pi, np.pi], np.array([1, 1]) / (2 * np.pi), label = "Uniform")
plt.title(file.replace("_", "\_").replace(".txt", ""))
plt.xlabel(r"$\phi$")
plt.ylabel(r"$P(\phi)$")
plt.legend()
plt.show()

# z coord is second
thetas = [np.arctan2(np.sqrt(_[0]**2 + _[1]**2), _[2]) for _ in data] 
plt.hist(thetas, density = True, bins = 100, range = (0, np.pi),
         label = "$N = " + str(len(thetas)) + "$")
xs = np.linspace(0, np.pi, 100)
# P(theta) is proportional to the circumfrence of a circle on sphere with given
# theta, so proportional to radius of that circle, so proportional to 
# sin(theta)
ys = np.sin(xs) / 2
plt.plot(xs, ys, label = r"$\sin(\theta)$")
plt.title(file.replace("_", "\_").replace(".txt", ""))
plt.xlabel(r"$\theta$")
plt.ylabel(r"$P(\theta)$")
plt.legend()
plt.show()

files = glob.glob("generate_point_on_plane_2d_*.txt")
files = sorted(files)
for file in files:
    data = np.loadtxt(file)
    h = float(file.split("_h_")[-1].split("_")[0])
    L = float(file.split("_L_")[-1][:-4])
    plt.hist(data[:, 0] / (L / 2.0), density = True, bins = 25, 
             range = (-1.0, 1.0), histtype = "step",
             label = "$h = " + str(h) + ", L = " + str(L) + "$")
    if not (data[:, 1] == h).all():
        print("Height failed!")
        exit()
    if (data[:, 0] < -L / 2.0).any():
        print("L failed!")
        exit()
    if (data[:, 0] >= L / 2.0).any():
        print("L failed!")
        exit()
plt.plot([-1.0, 1.0], np.array([0.5, 0.5]), label = "Uniform")
plt.title("generate point on plane 2d")
plt.xlabel(r"$x / (L / 2)$")
plt.ylabel(r"$P(x / (L / 2))$")
plt.legend(loc = "lower center", frameon = False)
plt.show()

files = glob.glob("generate_point_on_plane_3d_*.txt")
files = sorted(files)
for file in files:
    data = np.loadtxt(file)
    h = float(file.split("_h_")[-1].split("_")[0])
    L = float(file.split("_L_")[-1][:-4])
    if not (data[:, 1] == h).all():
        print("Height failed!")
        exit()
    symbols = ["x", "y", "z"]
    for i in [0, 2]:
        plt.hist(data[:, i] / (L / 2.0), density = True, bins = 25, 
                 range = (-1.0, 1.0), histtype = "step",
                 label = ("$h = " + str(h) + ", L = " + str(L) + ", " + 
                          symbols[i] + "$"))
        if (data[:, i] < -L / 2.0).any():
            print("L failed!")
            exit()
        if (data[:, i] >= L / 2.0).any():
            print("L failed!")
            exit()
plt.plot([-1.0, 1.0], np.array([0.5, 0.5]), label = "Uniform")
plt.title("generate point on plane 3d")
plt.xlabel(r"$x / (L / 2)$")
plt.ylabel(r"$P(x / (L / 2))$")
plt.legend(loc = "lower center", frameon = False)
plt.show()    

files = glob.glob("generate_jump_*.txt")
files = sorted(files)
for file in files:
    with open(file) as f:
        variance = float(f.readline().split("=")[-1])
    data = np.loadtxt(file, skiprows = 1)
    coordinates = ["X", "Y", "Z"]
    kDims = data.shape[1]
    var_sum = 0
    for i in range(kDims):
        var_sum += np.var(data[:, i])
        label = "$" + coordinates[i] + ", N = " + str(data.shape[0]) + "$"
        plt.hist(data[:, i], density = True, bins = 50, label = label, 
                 histtype = "step")
        jump_cutoff = float(file.split("_")[-2])
    print("jump_cutoff = {}".format(jump_cutoff))
    print("delta variance = {}".format(var_sum / kDims - variance))
    xs = np.linspace(-jump_cutoff, jump_cutoff, 100)
    dist = scipy.stats.norm()
    ys = dist.pdf(xs) / (dist.cdf(jump_cutoff) - dist.cdf(-jump_cutoff))
    plt.plot(xs, ys, label = r"pdf")
    plt.title(file.replace("_", "\_").replace(".txt", ""))
    plt.xlabel(r"$k$")
    plt.ylabel(r"$P(k)$")
    plt.legend()
    plt.show()

files = glob.glob("sample_first_hit_2d_*.txt")
files = sorted(files)
for i, file in enumerate(files):
    with open(file) as f:
        height = float(f.readline().split("=")[-1])
        particle = f.readline().split("=")[-1].split(",")[:-1]
        particle = np.array([float(_) for _ in particle])
    data = np.loadtxt(file)
    difference = particle[1] - height
    data = (data - particle) / difference
    X = 0
    Y = 1
    label = ("height = $" + str(height) + "$, particle = $" + 
             str(particle[X]) + ", " + str(particle[Y]) + "$")
    plt.hist(data[:, 0], density = True, bins = 100, range = (-10, 10),
                 histtype = "step",
                 label = label)
# plot analytical
xs = np.linspace(-10, 10, 1000)
ys =  1.0 / (2.0 * (xs**2 + 1.0)**(3.0/2.0))
plt.plot(xs, ys, label = "pdf", zorder = -1)
plt.xlabel(r"$x$")
plt.ylabel(r"$P(x)$")
plt.legend(loc = "upper left", frameon = False)
plt.title(r"sample first hit 2d")
plt.show()

files = glob.glob("sample_first_hit_3d_*.txt")
files = sorted(files)
for i, file in enumerate(files):
    with open(file) as f:
        height = float(f.readline().split("=")[-1])
        particle = f.readline().split("=")[-1].split(",")[:-1]
        particle = np.array([float(_) for _ in particle])
    data = np.loadtxt(file)
    difference = particle[1] - height
    data = (data - particle) / difference
    rs = np.sqrt(data[:, 0]**2 + data[:, 2]**2)
    X = 0
    Y = 1
    Z = 2
    label = ("height = $" + str(height) + "$, particle = $" + 
             str(particle[X]) + ", " + str(particle[Y]) + ", " + 
             str(particle[Z]) + "$")
    plt.figure(0)
    plt.hist(rs, density = True, bins = 100, range = (0, 10),
             histtype = "step", label = label)
    plt.figure(1)
    phis = np.arctan2(data[:, Z], data[:, X])
    plt.hist(phis, density = True, bins = 40, range = (-np.pi, np.pi),
                 histtype = "step", label = label)
# plot analytical
plt.figure(0)
xs = np.linspace(0.0, 10, 1000)
ys =  1.0 / ((xs**2 + 1.0)**(3.0/2.0))
plt.plot(xs, ys, zorder = -1, label = "pdf")
plt.xlabel(r"$r$")
plt.ylabel(r"$P(r)$")
plt.legend(loc = "upper right", frameon = False)
plt.title(r"sample first hit 3d")
plt.figure(1)
plt.plot([-np.pi, np.pi], np.array([1, 1]) / (2 * np.pi), label = "Uniform")
plt.xlabel(r"$\phi$")
plt.ylabel(r"$P(\phi)$")
plt.legend(loc = "lower center", frameon = False)
plt.title(r"sample first hit 3d")
plt.show()
