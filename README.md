# OpenPMDBeamphysics.jl

A Julia package for beam physics analysis following the
[openPMD-beamphysics](https://github.com/openPMD/openPMD-beamphysics) standard.
Adapted from the Python implementation by Christopher Mayes.

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/ndwang/openPMD-beamphysics.jl")
```

## Quick Start

```julia
using OpenPMDBeamphysics

pg = ParticleGroup("particles.h5")
pg["sigma_x"]       # weighted std of x [m]
pg["norm_emit_x"]   # normalized x emittance [m]
twiss(pg; plane="x")

fm = FieldMesh("solenoid.h5")
fm["Bz"]            # scaled Bz component
```

## Documentation

Full documentation including API reference is available at
https://ndwang.github.io/openPMD-beamphysics.
