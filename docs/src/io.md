# I/O

## File Format

This package reads and writes files in the
[openPMD-beamphysics](https://github.com/openPMD/openPMD-beamphysics) HDF5 format.

## Reading Particles

```julia
using OpenPMDBeamphysics

# Load particles (file must contain exactly one particle path)
pg = ParticleGroup("particles.h5")

# Inspect available paths in a file
using HDF5
h5open("particles.h5", "r") do h5
    println(particle_paths(h5))
end
```

See [`ParticleGroup`](@ref), [`particle_paths`](@ref).

## Writing Particles

```julia
write_particle_group("output.h5", pg)
```

See [`write_particle_group`](@ref).

## Reading Field Meshes

```julia
fm = FieldMesh("fieldmesh.h5")

# Inspect available paths
h5open("fieldmesh.h5", "r") do h5
    println(field_paths(h5))
end
```

See [`FieldMesh`](@ref), [`field_paths`](@ref).

## Writing Field Meshes

```julia
write_fieldmesh("output.h5", fm)
```

See [`write_fieldmesh`](@ref).

## Low-Level HDF5 API

For advanced use, the following functions are also exported. See the
[API Reference](@ref api) for full signatures.

**Reading:** [`particle_paths`](@ref), [`field_paths`](@ref),
[`load_field_attrs`](@ref), [`is_constant_component`](@ref),
[`constant_component_value`](@ref), [`component_data`](@ref),
[`offset_component_name`](@ref), [`particle_array`](@ref),
[`all_components`](@ref), [`component_str`](@ref)

**Writing:** [`pmd_init`](@ref), [`pmd_field_init`](@ref),
[`write_component_data`](@ref), [`write_pmd_bunch`](@ref),
[`write_pmd_field`](@ref), [`particle_data_dict`](@ref),
[`write_particle_group`](@ref)
