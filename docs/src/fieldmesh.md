# Field Mesh

## Overview

`FieldMesh` stores gridded electromagnetic field data in openPMD-beamphysics format.
It supports both rectangular (`x, y, z`) and cylindrical (`r, theta, z`) geometries,
static and RF fields, and complex field factors.

See [`FieldMesh`](@ref) in the [API Reference](@ref api) for full details.

## Loading and Saving

```julia
using OpenPMDBeamphysics

# Load from file (must contain exactly one field path)
fm = FieldMesh("solenoid.h5")

# Save to file
write_fieldmesh("output.h5", fm)
```

See [`FieldMesh`](@ref), [`write_fieldmesh`](@ref).

## Accessing Fields and Coordinates

Field components and coordinates are accessed via string indexing:

```julia
fm["Bz"]      # Bz component, scaled by field_scale * exp(i*phase)
fm["Er"]      # Er component (cylindrical only)
fm["z"]       # z coordinate vector
fm["r"]       # r coordinate vector (cylindrical)

# Operator prefixes (useful for complex RF fields)
fm["re_Bz"]   # real part of Bz
fm["im_Bz"]   # imaginary part of Bz
fm["abs_Ez"]  # |Ez|

# Field magnitudes
fm["B"]       # |B| across all components
fm["E"]       # |E| across all components
```

The underlying functions are [`coord_vec`](@ref) and [`scaled_component`](@ref);
use [`coord_vecs`](@ref) to retrieve all three coordinate vectors at once.

See [`coord_vec`](@ref), [`coord_vecs`](@ref), [`scaled_component`](@ref).

## Grid Information

```julia
mins(fm)            # grid minimum coordinates
maxs(fm)            # grid maximum coordinates
deltas(fm)          # grid spacing per axis
fm.grid_size        # number of grid points per axis
axis_index(fm, "z") # 1-based index of the "z" axis

# Per-axis helpers (rectangular)
xmin(fm); xmax(fm); dx(fm)
ymin(fm); ymax(fm); dy(fm)
zmin(fm); zmax(fm); dz(fm)

# Per-axis helpers (cylindrical)
rmin(fm); rmax(fm); dr(fm)
zmin(fm); zmax(fm); dz(fm)
```

See [`axis_index`](@ref), [`mins`](@ref), [`maxs`](@ref), [`deltas`](@ref).

## Field Properties

```julia
is_static(fm)        # true if harmonic == 0
is_pure_magnetic(fm) # true if no electric components
is_pure_electric(fm) # true if no magnetic components
scale(fm)            # field_scale factor
phase(fm)            # phase in radians: -2π * RFphase
factor(fm)           # complex factor: scale * exp(i*phase)
frequency(fm)        # oscillation frequency [Hz]
```

See [`is_static`](@ref), [`is_pure_magnetic`](@ref), [`is_pure_electric`](@ref),
[`scale`](@ref), [`phase`](@ref), [`factor`](@ref), [`frequency`](@ref).

## Grid Modification

```julia
set_min!(fm, "z", -0.1)  # shift z origin to -0.1 m
set_max!(fm, "z",  0.1)  # shift z origin so max = 0.1 m
```

See [`set_min!`](@ref), [`set_max!`](@ref).
