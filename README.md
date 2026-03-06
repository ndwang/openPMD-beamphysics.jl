# BeamPhysics.jl

A Julia package for beam physics analysis following the
[openPMD-beamphysics](https://github.com/openPMD/openPMD-standard/blob/upcoming-2.0.0/EXT_BeamPhysics.md) standard.
Inspired by the Python package [openPMD-beamphysics](https://github.com/ChristopherMayes/openPMD-beamphysics).

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/ndwang/BeamPhysics.jl")
```

## Quick Start

```julia
using BeamPhysics

pg = ParticleGroup("particles.h5")
pg["sigma_x"]       # weighted std of x [m]
pg["norm_emit_x"]   # normalized x emittance [m]
twiss(pg; plane="x")

fm = FieldMesh("solenoid.h5")
fm["Bz"]            # scaled Bz component
```

## Documentation

Full documentation including API reference is available at
https://ndwang.github.io/BeamPhysics.jl.
