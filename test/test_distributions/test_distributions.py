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

files = glob.glob("sample_first_hit_*_2d.txt")
files = sorted(files)
alphas = []
for i, file in enumerate(files):
    with open(file) as f:
        line = f.readline()
        particle_x = float(line.split()[-2])
        particle_y = float(line.split()[-1])
    particle = np.array([particle_x, particle_y])
    alpha = np.linalg.norm(particle)
    alphas.append(alpha)
    X = 0
    Y = 1
    particle_angle = np.arctan2(particle[Y], particle[X])
    label = "particle = " + str(particle[X]) + ", " + str(particle[Y])
    data = np.loadtxt(file, skiprows = 1)
    # y coord is first
    angles = [np.arctan2(_[1], _[0]) - particle_angle for _ in data] 
    angles = np.array(angles)
    angles[angles < -np.pi] += 2 * np.pi
    plt.hist(angles, density = True, bins = 50, 
             range = (-np.pi, np.pi), label = label, histtype = "step")
alphas = set(alphas)
for alpha in alphas:
    xs = np.linspace(-np.pi, np.pi, 100)
    ys = (alpha**2 - 1) / (alpha**2 - 2 * alpha * np.cos(xs) + 1) 
    ys = ys / (2 * np.pi)
    plt.plot(xs, ys, label = r"pdf, $\alpha = " + "{:.2f}".format(alpha) + "$")
plt.xlabel(r"$\theta$")
plt.ylabel(r"$P(\theta)$")
plt.legend(loc = "upper left")
plt.title(r"sample\_first\_hit 2d, $N = " + str(len(angles)) + "$")
plt.show()

# projects v into plane containing phi_vector perpendicular to particle 
# and gets angle (phi) relative to phi_vector
def get_phi(particle, phi_vector, v):
    particle = particle / np.linalg.norm(particle)
    v = v - np.dot(v, particle) * particle
    v = v / np.linalg.norm(v)
    x = np.dot(v, phi_vector)
    # taking dot product with y axis defined by cross product
    y = np.dot(v - np.dot(v, phi_vector) * phi_vector, 
               np.cross(particle, phi_vector))
    phi = np.arctan2(y, x)
    return phi

files = glob.glob("sample_first_hit_*_3d.txt")
files = sorted(files)
alphas = []
for i, file in enumerate(files):
    with open(file) as f:
        line = f.readline()
        particle_x = float(line.split()[-3])
        particle_y = float(line.split()[-2])
        particle_z = float(line.split()[-1])
    particle = np.array([particle_x, particle_y, particle_z])
    alpha = np.linalg.norm(particle)
    alphas.append(alpha)
    X = 0
    Y = 1
    Z = 2
    # need to define phi, generate vector perpendicular to parti
    tmp = np.array([0, 0, 1])
    if np.abs(np.dot(tmp, particle / np.linalg.norm(particle))) < 1e-3:
        raise Exception("particle cannot point along z axis for these tests!")
    data = np.loadtxt(file, skiprows = 1)
    data = [_ / np.linalg.norm(_) for _ in data]
    # define phi relative to this vector
    phi_vector = np.cross(particle, tmp)
    phi_vector = phi_vector / np.linalg.norm(phi_vector)
    phis = [get_phi(particle, phi_vector, _) for _ in data]
    label = ("particle = " + str(particle[X]) + ", " + str(particle[Y]) + 
             ", " + str(particle[Z]))
    plt.figure(1)
    plt.hist(phis, density = True, bins = 50, range = (-np.pi, np.pi), 
             label = label, histtype = "step")
    plt.figure(2)
    cos_thetas = [np.dot(particle / np.linalg.norm(particle), _) for _ in data]
    cos_thetas = np.array(cos_thetas)
    plt.hist(cos_thetas, density = True, bins = 50, range = (-1, 1), 
             label = label, histtype = "step")
plt.figure(1)
plt.plot([-np.pi, np.pi], np.array([1, 1]) / (2 * np.pi), label = "Uniform")
plt.xlabel(r"$\phi$")
plt.ylabel(r"$P(\phi)$")
plt.legend(loc = "best")
plt.title(r"sample\_first\_hit 3d, $N = " + str(len(phis)) + "$")
plt.figure(2)
alphas = set(alphas)
for alpha in alphas:
    # analytical pdf
    xs = np.linspace(-1, 1, 100)
    ys = (alpha / 2 * (alpha**2 - 1) / 
          (alpha**2 - 2 * alpha * xs + 1)**(3.0 / 2.0))
    # need to add uniform component back in if particle escapes
    ys = 1 / alpha * ys + (1 - 1 / alpha) * (1.0 / 2.0)
    plt.plot(xs, ys, label = r"pdf, $\alpha = " + "{:.2f}".format(alpha) + "$")
#plt.plot([-np.pi, np.pi], np.array([1, 1]) / (2 * np.pi), label = "Uniform")
plt.xlabel(r"$\cos(\theta)$")
plt.ylabel(r"$P(\cos(\theta))$")
plt.legend(loc = "best")
plt.title(r"sample\_first\_hit 3d, $N = " + str(len(phis)) + "$")
plt.show()

    
