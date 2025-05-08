"""
Utilities for writing particle and field data to HDF5 files in openPMD format.
"""

"""
    pmd_init(h5, basePath="/data/%T/", particlesPath="./")

Root attribute initialization.

h5 should be the root of the file.
"""
function pmd_init(h5, basePath="/data/%T/", particlesPath="./")
    d = Dict(
        "basePath" => basePath,
        "dataType" => "openPMD",
        "openPMD" => "2.0.0",
        "openPMDextension" => "BeamPhysics;SpeciesType",
        "particlesPath" => particlesPath
    )
    
    for (k, v) in d
        h5.attrs[k] = fstr(v)
    end
end

"""
    pmd_field_init(h5, externalFieldPath="/ExternalFieldPath/%T/")

Root attribute initialization for an openPMD-beamphysics External Field Mesh

h5 should be the root of the file.
"""
function pmd_field_init(h5, externalFieldPath="/ExternalFieldPath/%T/")
    d = Dict(
        "dataType" => "openPMD",
        "openPMD" => "2.0.0",
        "openPMDextension" => "BeamPhysics",
        "externalFieldPath" => externalFieldPath
    )
    
    for (k, v) in d
        h5.attrs[k] = fstr(v)
    end
end

"""
    write_pmd_bunch(h5, data, name=nothing)

Data is a dict with:
    Array: 'x', 'px', 'y', 'py', 'z', 'pz', 't', 'status', 'weight'
    String: 'species'
    Int: n_particle

Optional data:
    Array: 'id'

See inverse routine:
    .particles.load_bunch_data
"""
function write_pmd_bunch(h5, data, name=nothing)
    g = isnothing(name) ? h5 : create_group(h5, name)

    # Write into species group
    species = data["species"]
    g = create_group(g, species)

    # Attributes
    g.attrs["speciesType"] = fstr(species)
    g.attrs["numParticles"] = data["n_particle"]
    g.attrs["totalCharge"] = data["charge"]
    g.attrs["chargeUnitSI"] = 1.0

    # Required Datasets
    for key in ["x", "px", "y", "py", "z", "pz", "t", "status", "weight"]
        # Get full name, write data
        g2_name = component_from_alias[key]

        # Units
        u = pg_units(key)

        # Write
        write_component_data(g, g2_name, data[key], unit=u)
    end

    # Optional id. This does not have any units.
    if haskey(data, "id")
        g["id"] = data["id"]
    end
end

"""
    write_pmd_field(h5, data, name=nothing)

Data is a dict with:
    attrs: flat dict of attributes.
    components: flat dict of components

See inverse routine:
    .readers.load_field_data
"""
function write_pmd_field(h5, data, name=nothing)
    g = isnothing(name) ? h5 : create_group(h5, name)

    # Validate attrs
    attrs, other = load_field_attrs(data["attrs"])

    # Encode for writing
    attrs = encode_attrs(attrs)

    # Write attributes
    for (k, v) in attrs
        g.attrs[k] = v
    end

    # All other attributes (don't change them)
    for (k, v) in other
        g.attrs[k] = v
    end

    # write components (datasets)
    for (key, val) in data["components"]
        # Units
        u = pg_units(key)

        # Ensure complex
        val = complex.(val)

        # Write
        write_component_data(g, key, val, unit=u)
    end
end

"""
    write_component_data(h5, name, data, unit=nothing)

Writes data to a dataset h5[name]

If data is a constant array, a group is created with the constant value and shape

If unit is given, this will be used
"""
function write_component_data(h5, name, data, unit=nothing)
    # Check for constant component
    dat0 = data[1]

    if all(x -> x == dat0, data)
        g = create_group(h5, name)
        g.attrs["value"] = dat0
        g.attrs["shape"] = size(data)
    else
        h5[name] = data
        g = h5[name]
    end

    if !isnothing(unit)
        g.attrs["unitSI"] = unit.unitSI
        g.attrs["unitDimension"] = unit.unitDimension
        g.attrs["unitSymbol"] = fstr(unit.unitSymbol)
    end

    return g
end