# Particles

## Overview

`ParticleGroup` is the primary data structure for beam particles. It stores particle
coordinates, momenta, time, status, weight, and species. Array fields are parameterized
so they can be `Vector`, `CuArray`, or any `AbstractVector`, enabling GPU computation.

See [`ParticleGroup`](@ref) in the [API Reference](@ref api) for the full type definition.

## Creating a ParticleGroup

```julia
using OpenPMDBeamphysics

# Load from an openPMD HDF5 file
pg = ParticleGroup("particles.h5")

# Create a single test particle
pg = single_particle(z=0.01, pz=1e6)

# Construct directly
using AtomicAndPhysicalConstants: Species
pg = ParticleGroup(x, px, y, py, z, pz, t, status, weight, Species("electron"))
```

See [`single_particle`](@ref), [`ParticleGroup`](@ref) (file or HDF5 group constructors).

## Accessing Properties

Stored fields and derived properties can be accessed as **direct function calls**
or via **string-based indexing** — both are equivalent:

```julia
# Direct function calls
energy(pg)       # total energy array [eV]
gamma(pg)        # Lorentz factor array
r(pg)            # cylindrical radius array [m]
norm_emit_x(pg)  # normalized x emittance scalar [m]
charge(pg)       # total charge scalar [C]

# Equivalent string-based indexing
pg["energy"]
pg["gamma"]
pg["r"]
pg["norm_emit_x"]
pg["charge"]
```

Stored fields are accessible as struct fields (`pg.x`, `pg.pz`, etc.) or via `pg["x"]`.

String indexing additionally supports statistical prefixes that have no direct
function equivalent:

```julia
pg["sigma_x"]    # std(x), weighted
pg["mean_pz"]    # mean(pz), weighted
pg["min_x"]      # minimum x
pg["max_x"]      # maximum x
pg["ptp_z"]      # peak-to-peak range of z
pg["cov_x_px"]   # covariance of x and px
pg["delta_x"]    # x .- mean(x) per particle
```

## Slicing

```julia
pg[1:100]        # first 100 particles (returns a copy)
view(pg, 1:100)  # non-owning view into the first 100 particles
```

See [`ParticleGroupView`](@ref).

## Joining

```julia
pg3 = pg1 + pg2  # concatenate two ParticleGroups of the same species
```

## In-place Operations

```julia
drift!(pg, 1e-9)       # drift by 1 ns
drift_to_z!(pg, 0.5)   # drift all particles to z = 0.5 m
drift_to_t!(pg)        # drift all particles to their mean time
set_charge!(pg, 1e-12) # rescale weights to total charge = 1 pC
```

## Splitting

```julia
chunks = split_particles(pg; n_chunks=100, key="z")
```

## Derived Properties

All derived properties are callable as `f(pg)` **and** accessible via `pg["key"]`.

| Function call | String key | Returns | Units |
|---|---|---|---|
| `nalive(pg)` | `pg["n_alive"]` | number of alive particles | — |
| `ndead(pg)` | `pg["n_dead"]` | number of dead particles | — |
| `charge(pg)` | `pg["charge"]` | total charge | C |
| `mass(pg)` | `pg["mass"]` | particle rest mass | eV |
| `species_charge(pg)` | `pg["species_charge"]` | elementary charge of species | C |
| `average_current(pg)` | `pg["average_current"]` | average current | A |
| `momentum(pg)` | `pg["p"]` | total momentum per particle | eV/c |
| `energy(pg)` | `pg["energy"]` | total energy per particle | eV |
| `kinetic_energy(pg)` | `pg["kinetic_energy"]` | kinetic energy per particle | eV |
| `gamma(pg)` | `pg["gamma"]` | Lorentz factor per particle | — |
| `beta(pg)` | `pg["beta"]` | `\|v\|/c` per particle | — |
| `beta_x(pg)` | `pg["beta_x"]` | `vx/c` per particle | — |
| `beta_y(pg)` | `pg["beta_y"]` | `vy/c` per particle | — |
| `beta_z(pg)` | `pg["beta_z"]` | `vz/c` per particle | — |
| `xp(pg)` | `pg["xp"]` | `px/pz` per particle | rad |
| `yp(pg)` | `pg["yp"]` | `py/pz` per particle | rad |
| `r(pg)` | `pg["r"]` | cylindrical radius per particle | m |
| `theta(pg)` | `pg["theta"]` | azimuthal angle per particle | rad |
| `pr(pg)` | `pg["pr"]` | radial momentum per particle | eV/c |
| `ptheta(pg)` | `pg["ptheta"]` | azimuthal momentum per particle | eV/c |
| `Lz(pg)` | `pg["Lz"]` | angular momentum per particle | m·eV/c |
| `norm_emit_x(pg)` | `pg["norm_emit_x"]` | normalized x emittance | m |
| `norm_emit_y(pg)` | `pg["norm_emit_y"]` | normalized y emittance | m |
| `norm_emit_4d(pg)` | `pg["norm_emit_4d"]` | normalized 4D emittance | m² |
| `x_bar(pg)` | `pg["x_bar"]` | normalized x coordinate per particle | m^(1/2) |
| `px_bar(pg)` | `pg["px_bar"]` | normalized px coordinate per particle | m^(1/2) |
| `y_bar(pg)` | `pg["y_bar"]` | normalized y coordinate per particle | m^(1/2) |
| `py_bar(pg)` | `pg["py_bar"]` | normalized py coordinate per particle | m^(1/2) |
| `Jx(pg)` | `pg["Jx"]` | x action amplitude per particle | m |
| `Jy(pg)` | `pg["Jy"]` | y action amplitude per particle | m |
