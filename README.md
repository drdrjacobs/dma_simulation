# Reactive Deposition Sim

A simulation of Brownian particles depositing from all sides onto a uniformly reactive cluster. The code can be run in both two and three dimensions. 

Initially, the cluster displays a compact circular/spherical structure. However, upon reaching a critical radius, it forms dendritic branches. This transition is illustrated in the figure below.

This kind of dendritic deposition occurs on the anode in next-generation lithium-metal batteries, and the problems it creates represent one of the main barriers to commercializing the technology. The fractal pattern formed by the cluster also appears in many other contexts, including lightning strikes, Hele-Shaw viscous fingering, and even bacterial growth! 

For more information, see:

Jacobson, D. & Miller III, T. F. Compact-to-Dendritic Transition in the
Reactive Deposition of Brownian Particles. [arXiv:2204.01173 [cond-mat]](https://arxiv.org/abs/2204.01173) (2022).

| <img src="https://raw.githubusercontent.com/drdrjacobs/reactive_deposition_sim/master/images/dendritic_transition.png" width="600" height="803"> | 
|:--:| 
| The compact-to-dendritic transition in the reactive deposition of Brownian particles at $\log_{10} \mathrm{Da} = -2.12$. Color schemes vary depending on the panel, as indicated in the descriptions. Panels (a) and (b) depict two-dimensional deposition. The $N$ particles in these panels are rendered with a radius twice that of the actual radius for clarity. The $i$-th band of color moving outwards from the center corresponds to the structure after $i N / 5$ particles have been deposited. (a) Initially, the cluster displays a compact morphology. (b) However, upon reaching a critical radius, it spontaneously forms dendritic branches. (c) Main Panel: During three-dimensional deposition, a similar dendritic morphology emerges at large length scales. Particles are colored based on their distance from the origin. Inset: A two-dimensional slice through the initial compact 3D cluster, $N = 1 \times 10^{6}$. The color scheme follows a continuous gradient based on the deposition order. |

# Dependencies

- [Boost](https://www.boost.org/)
- [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page): Template library for linear algebra.

# Usage

### Config File

First generate a config file for Make starting from `example_config.mk`:

```cp example_config.mk config.mk```

Change `BOOST_PATH` and `EIGEN_INC_PATH` in `config.mk` to point to your installations of Boost and Eigen.

### Compile

Create the build directories:

```mkdir obj2d obj3d bin```

And compile with:

```make```

This will generete two executables, `bin/ga_simulation_2d` for running two-dimensional simulations and `bin/ga_simulation_3d` for running three-dimensional simulations.

### Run

Finally, modify the `params.txt` file to specify physical system you would like to simulate. You can now run with:

```./bin/ga_simulation_2d```

After, the simulation finishes, it will output a `frame_NNN.xyz` file (where `NNN` is the number of particles) that lists the positions off all the particles in the cluster. This file can be vizualized with [Ovito](https://www.ovito.org/).


