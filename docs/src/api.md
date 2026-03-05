# [API Reference](@id api)

Complete listing of all exported symbols.

## Types

```@docs
ParticleGroup
ParticleGroupView
FieldMesh
PMDUnit
```

## Particle Constructors

```@docs
single_particle
ParticleGroup(h5::HDF5.Group)
ParticleGroup(file::String)
```

## Particle Derived Properties

The following functions are available both as standalone calls and via string
indexing (`pg["key"]`). Most are undocumented one-liners; see source for details.

**Scalar:** `nalive`, `ndead`, `mass`, `species_charge`, `charge`, `average_current`

**Per-particle arrays:**
`momentum`, `energy`, `kinetic_energy`, `gamma`, `beta`, `beta_x`, `beta_y`, `beta_z`,
`xp`, `yp`, `r`, `theta`, `pr`, `ptheta`, `Lz`,
`x_bar`, `px_bar`, `y_bar`, `py_bar`, `Jx`, `Jy`,
`norm_emit_x`, `norm_emit_y`, `norm_emit_4d`

## Particle Operations

`drift!(pg, dt)`, `drift_to_z!(pg[, z])`, `drift_to_t!(pg[, t])`,
`set_charge!(pg, val)`, `split_particles(pg; n_chunks, key)`,
`in_z_coordinates(pg)`, `in_t_coordinates(pg)`

## Statistics

```@docs
norm_emit_calc
twiss_calc
twiss_dispersion
twiss
A_mat_calc
A_inverse_mat_calc
amplitude_calc
particle_amplitude
normalized_particle_coordinate
slice_statistics
```

## FieldMesh Accessors

```@docs
axis_index
mins
maxs
deltas
coord_vec
coord_vecs
dx
dy
dz
dr
dtheta
xmin
xmax
ymin
ymax
zmin
zmax
rmin
rmax
thetamin
thetamax
scale
phase
factor
frequency
is_static
is_pure_electric
is_pure_magnetic
scaled_component
set_min!
set_max!
write_fieldmesh
```

## Units & Display

```@docs
PMDUnit(symbol::String)
nice_scale_prefix
nice_array
limits
plottable_array
pg_units
mathlabel
texlabel
```

## I/O

```@docs
particle_paths
field_paths
load_field_attrs
is_constant_component
constant_component_value
component_data
offset_component_name
particle_array
all_components
component_str
pmd_init
pmd_field_init
write_component_data
write_pmd_bunch
write_pmd_field
particle_data_dict
write_particle_group
```

## Plotting

```@docs
density_plot
marginal_plot
slice_plot
density_and_slice_plot
fieldmesh_plot
fieldmesh_plot_onaxis
```
