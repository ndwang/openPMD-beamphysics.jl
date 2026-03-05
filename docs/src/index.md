# BeamPhysics.jl

A Julia package for working with beam physics data in the
[openPMD-beamphysics](https://github.com/openPMD/openPMD-beamphysics) format.

## Features

- Load and save particle data (`ParticleGroup`) from openPMD HDF5 files
- Load and save gridded electromagnetic field data (`FieldMesh`) from openPMD HDF5 files
- Compute derived particle properties (emittance, Twiss parameters, etc.)
- Drift particles in free space
- Optional plotting via [Plots.jl](https://github.com/JuliaPlots/Plots.jl)

## Installation

```julia
using Pkg
Pkg.add("BeamPhysics")
```

## Quick Start

```julia
using BeamPhysics

# Load particles from an openPMD HDF5 file
pg = ParticleGroup("particles.h5")

# Access particle properties
pg["x"]           # x positions [m]
pg["sigma_x"]     # std of x [m]
pg["norm_emit_x"] # normalized emittance [m]

# Load a field mesh
fm = FieldMesh("solenoid.h5")
fm["Bz"]          # scaled Bz component
coord_vec(fm, "z") # z coordinate vector
```

## Contents

```@contents
Pages = ["particles.md", "fieldmesh.md", "statistics.md", "io.md", "plotting.md", "api.md"]
Depth = 2
```
