# BeamTracking Integration

BeamPhysics.jl provides a package extension for
[BeamTracking.jl](https://github.com/bmad-sim/BeamTracking.jl) that enables
direct conversion from `ParticleGroup` to `Bunch`.

The extension is loaded automatically when both packages are imported.

## Usage

```julia
using BeamPhysics
using BeamTracking

# From a ParticleGroup
pg = ParticleGroup("particles.h5")
bunch = Bunch(pg)

# Directly from an openPMD HDF5 file
bunch = Bunch("particles.h5")
```

## Coordinate Mapping

The `ParticleGroup` stores absolute coordinates (positions in meters, momenta in
eV/c, time in seconds). The `Bunch` uses Bmad-style normalized phase space
coordinates relative to a reference particle. The conversion maps:

| Bunch coordinate | Formula | Description |
|---|---|---|
| `v[:,1]` (x) | `pg.x` | Horizontal position [m] |
| `v[:,2]` (px) | `pg.px / p_ref` | Normalized horizontal momentum |
| `v[:,3]` (y) | `pg.y` | Vertical position [m] |
| `v[:,4]` (py) | `pg.py / p_ref` | Normalized vertical momentum |
| `v[:,5]` (z) | `(pg.z - z_ref) - β_ref c (pg.t - t_ref)` | Longitudinal offset [m] |
| `v[:,6]` (pz) | `\|p\| / p_ref - 1` | Fractional momentum deviation δ |

where `p_ref` is the reference total momentum (mean of alive particles),
`z_ref` is the mean z of alive particles, `t_ref` is the mean time of alive
particles, and `β_ref = p_ref / √(p_ref² + (mc²)²)`.

The `z` formula reduces to:
- `z - z_ref` when all particles share the same time (t-coordinates)
- `-β_ref c (t - t_ref)` when all particles share the same z (z-coordinates)

## Reference Values

By default, `p_over_q_ref` and `t_ref` are computed from the alive particles.
They can be overridden:

```julia
bunch = Bunch(pg; p_over_q_ref=my_R, t_ref=0.0)
```

`p_over_q_ref` follows the BeamTracking convention: `p_ref / (c * q)` where `q`
is the species charge from `chargeof(species)`.

## Particle Status

Particle status values are copied directly. `ParticleGroup` uses `1` for alive
particles, which matches BeamTracking's `STATE_ALIVE`.

## Limitations

- Spin information is not converted (`ParticleGroup` has no spin fields); the
  `Bunch` quaternion field is set to `nothing`.
- Particle weights (macro-charge) from `ParticleGroup` are copied to the `Bunch`
  `weight` field.
