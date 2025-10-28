# openPMD-beamphysics

A Julia package for beam physics analysis using the openPMD standard.

## Overview

openPMD-beamphysics is a Julia package for beam physics analysis and visualization. The package follows the [openPMD standard BeamPhysics extension](https://github.com/openPMD/openPMD-standard/blob/upcoming-2.0.0/EXT_BeamPhysics.md) for data storage and exchange, making it compatible with various particle accelerator simulation codes. It provides tools for particle tracking, statistical analysis, emittance calculations, and visualization of beam dynamics.

This package is adapted from the Python implementation [openPMD-beamphysics](https://github.com/ChristopherMayes/openPMD-beamphysics) and provides native Julia performance with full compatibility with the openPMD data format.

## Installation

To install the package
```julia
using Pkg
Pkg.add(url="https://github.com/ndwang/openPMD-beamphysics.jl")
```

## Features

### Core Functionality
- **Particle Group Management**: Handle collections of particles with position, momentum, and metadata
- **Statistical Analysis**: Calculate beam parameters including emittance, Twiss parameters, and beam moments
- **openPMD Compatibility**: Read and write data in the openPMD standard format
- **Visualization**: Comprehensive plotting tools for beam analysis

## Quick Start

```julia
using openPMD_beamphysics

# Create a simple particle group
pg = single_particle(x=1e-3, px=1e-6, y=0.0, py=0.0, z=0.0, pz=0.0)

# Calculate emittance
emit_x = norm_emit_calc(pg, ["x"])

# Calculate Twiss parameters
sigma_mat = cov(pg, "x", "px")
twiss_params = twiss_calc(sigma_mat)

# Plot beam properties
slice_plot(pg, "sigma_x", "sigma_px")
```

## Main Types and Functions

### Core Types

- `ParticleGroup`: Main type representing a collection of particles

### Key Functions

#### Statistical Analysis
- `norm_emit_calc(particle_group, planes)`: Calculate normalized emittance
- `twiss_calc(sigma_mat)`: Calculate Twiss parameters from covariance matrix
- `twiss_dispersion(particle_group)`: Calculate dispersion parameters
- `slice_statistics(particle_group, keys...)`: Calculate statistics in slices

#### Particle Operations
- `single_particle(;kwargs...)`: Create a single particle ParticleGroup
- `drift!(particle_group, distance)`: Drift particles
- `split_particles(particle_group, n)`: Split particles for analysis

#### Visualization
- `slice_plot(particle_group, keys...)`: Create slice plots
- `density_plot(particle_group, x_key, y_key)`: Create density plots
- `marginal_plot(particle_group, keys...)`: Create marginal distribution plots

## Examples

### Basic Beam Analysis

```julia
using openPMD_beamphysics

# Create a Gaussian beam
n_particles = 10000
x = randn(n_particles) * 1e-3  # 1 mm RMS
px = randn(n_particles) * 1e-6  # 1 μrad RMS
y = randn(n_particles) * 1e-3
py = randn(n_particles) * 1e-6
z = randn(n_particles) * 1e-2  # 1 cm RMS
pz = randn(n_particles) * 1e-3

pg = ParticleGroup(x, px, y, py, z, pz, zeros(n_particles), 
                  ones(Int, n_particles), ones(n_particles), 
                  Species("electron"))

# Calculate beam parameters
emit_x = norm_emit_calc(pg, ["x"])
emit_y = norm_emit_calc(pg, ["y"])
emit_4d = norm_emit_calc(pg, ["x", "y"])

println("X emittance: $(emit_x) m⋅rad")
println("Y emittance: $(emit_y) m⋅rad")
println("4D emittance: $(emit_4d) m²⋅rad²")
```

### Twiss Parameter Analysis

```julia
# Calculate covariance matrix
sigma_mat = cov(pg, "x", "px")

# Get Twiss parameters
twiss = twiss_calc(sigma_mat)
println("Beta: $(twiss["beta"]) m")
println("Alpha: $(twiss["alpha"])")
println("Gamma: $(twiss["gamma"]) 1/m")
```

### Visualization

```julia

# Create slice plot
p1 = slice_plot(pg, "sigma_x", "sigma_px", n_slice=20)

# Create density plot
p2 = density_plot(pg, "x", "px")

# Combine plots
plot(p1, p2, layout=(1,2), size=(800,400))
```


## Authors

- ndwang <ndwang3@gmail.com>
