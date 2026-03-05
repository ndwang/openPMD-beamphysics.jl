# Plotting

Plotting functionality is available as a package extension when
[Plots.jl](https://github.com/JuliaPlots/Plots.jl) is loaded.

## Setup

```julia
using OpenPMDBeamphysics
using Plots
```

## Particle Plots

```julia
# 2D density plot of two particle properties
density_plot(pg, "x", "y")
density_plot(pg, "z", "energy"; nbins=100)

# Marginal plot: 2D density with 1D projections on margins
marginal_plot(pg, "z", "x")

# Slice plot: statistics along z slices
slice_plot(pg, "z", "sigma_x")

# Combined density + slice
density_and_slice_plot(pg, "z", "x")
```

See [`density_plot`](@ref), [`marginal_plot`](@ref), [`slice_plot`](@ref),
[`density_and_slice_plot`](@ref).

## Field Mesh Plots

```julia
# 2D field plot (slice through the mesh)
fieldmesh_plot(fm, "Bz")
fieldmesh_plot(fm, "Bz"; ix=1)  # specify slice index

# On-axis field profile
fieldmesh_plot_onaxis(fm, "Bz")
fieldmesh_plot_onaxis(fm, "Ez"; ir=1)
```

See [`fieldmesh_plot`](@ref), [`fieldmesh_plot_onaxis`](@ref).
