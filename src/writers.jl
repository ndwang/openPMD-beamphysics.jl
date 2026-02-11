"""
Utilities for writing particle and field data to HDF5 files in openPMD format.
"""

"""
    pmd_init(h5; basePath="/data/%T/", particlesPath="./")

Root attribute initialization for openPMD particle files.
h5 should be the root of the file.
"""
function pmd_init(h5; basePath="/data/%T/", particlesPath="./")
    attrs(h5)["basePath"] = basePath
    attrs(h5)["dataType"] = "openPMD"
    attrs(h5)["openPMD"] = "2.0.0"
    attrs(h5)["openPMDextension"] = "BeamPhysics;SpeciesType"
    attrs(h5)["particlesPath"] = particlesPath
end

"""
    pmd_field_init(h5; externalFieldPath="/ExternalFieldPath/%T/")

Root attribute initialization for an openPMD-beamphysics External Field Mesh.
h5 should be the root of the file.
"""
function pmd_field_init(h5; externalFieldPath="/ExternalFieldPath/%T/")
    attrs(h5)["dataType"] = "openPMD"
    attrs(h5)["openPMD"] = "2.0.0"
    attrs(h5)["openPMDextension"] = "BeamPhysics"
    attrs(h5)["externalFieldPath"] = externalFieldPath
end

"""
    write_component_data(h5, name, data; unit=nothing)

Write a single record component to an HDF5 group.

If all values in `data` are equal (constant component), a group is created with
`value` and `shape` stored as attributes. Otherwise, a dataset is created.

Unit metadata is written via `write_unit_h5` when `unit` is provided.
"""
function write_component_data(h5, name, data; unit=nothing)
    dat0 = data[1]

    if all(x -> x == dat0, data)
        # Constant component: store as group with value/shape attributes
        g = create_group(h5, name)
        attrs(g)["value"] = dat0
        attrs(g)["shape"] = collect(size(data))
    else
        # Non-constant: store as dataset
        h5[name] = data
        g = h5[name]
    end

    if !isnothing(unit)
        write_unit_h5(g, unit)
    end

    return g
end

"""
    write_pmd_bunch(h5, data::Dict; name=nothing)

Write particle bunch data in openPMD format.

`data` is a Dict with:
- Array keys: `"x"`, `"px"`, `"y"`, `"py"`, `"z"`, `"pz"`, `"t"`, `"status"`, `"weight"`
- String key: `"species"`
- Int key: `"n_particle"`
- Float key: `"charge"`
- Optional array key: `"id"`

See inverse routine: `ParticleGroup(h5::HDF5.Group)`
"""
function write_pmd_bunch(h5, data::Dict; name=nothing)
    g = isnothing(name) ? h5 : create_group(h5, name)

    # Create species subgroup
    species = data["species"]
    g = create_group(g, species)

    # Species attributes
    attrs(g)["speciesType"] = species
    attrs(g)["numParticles"] = data["n_particle"]
    attrs(g)["totalCharge"] = data["charge"]
    attrs(g)["chargeUnitSI"] = 1.0

    # Required components
    for key in ["x", "px", "y", "py", "z", "pz", "t", "status", "weight"]
        g2_name = COMPONENT_FROM_ALIAS[key]
        u = pg_units(key)
        write_component_data(g, g2_name, data[key]; unit=u)
    end

    # Optional id (no units)
    if haskey(data, "id")
        g["id"] = data["id"]
    end
end

"""
    write_pmd_field(h5, data::Dict; name=nothing)

Write field mesh data in openPMD format.

`data` is a Dict with:
- `"attrs"`: flat dict of field attributes
- `"components"`: flat dict of component name => array

See inverse routine: `load_field_attrs`
"""
function write_pmd_field(h5, data::Dict; name=nothing)
    g = isnothing(name) ? h5 : create_group(h5, name)

    # Validate attrs
    field_attrs, other = load_field_attrs(data["attrs"])

    # Write attributes
    for (k, v) in field_attrs
        attrs(g)[k] = v
    end
    for (k, v) in other
        attrs(g)[k] = v
    end

    # Write components
    for (key, val) in data["components"]
        u = pg_units(key)
        val = complex.(val)
        write_component_data(g, key, val; unit=u)
    end
end

"""
    particle_data_dict(pg::ParticleGroup) -> Dict{String, Any}

Extract a Dict from a ParticleGroup in the format expected by `write_pmd_bunch`.
"""
function particle_data_dict(pg::ParticleGroup)
    return Dict{String, Any}(
        "x" => pg.x,
        "px" => pg.px,
        "y" => pg.y,
        "py" => pg.py,
        "z" => pg.z,
        "pz" => pg.pz,
        "t" => pg.t,
        "status" => pg.status,
        "weight" => pg.weight,
        "species" => nameof(pg.species),
        "n_particle" => length(pg),
        "charge" => charge(pg),
        "id" => pg.id,
    )
end

"""
    write_particle_group(filename::String, pg::ParticleGroup)

Write a ParticleGroup to an HDF5 file in openPMD format.
This is the inverse of `ParticleGroup(filename::String)`.
"""
function write_particle_group(filename::String, pg::ParticleGroup)
    h5open(filename, "w") do h5
        pmd_init(h5)
        g = create_group(h5, "data/1")
        write_pmd_bunch(g, particle_data_dict(pg))
    end
end
