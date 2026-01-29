"""
Utilities for writing particle and field data to HDF5 files in openPMD format.
"""

# String formatting for HDF5 attributes (identity for now)
fstr(s::String) = s

# Encode field attributes for HDF5 writing (identity for now)
encode_attrs(attrs::Dict) = attrs

"""
    pmd_init(h5, basePath="/data/%T/", particlesPath="./")

Root attribute initialization.

h5 should be the root of the file.
"""
function pmd_init(h5, basePath="/data/%T/", particlesPath="./")
    h5["basePath"] = fstr(basePath)
    h5["dataType"] = fstr("openPMD")
    h5["openPMD"] = fstr("2.0.0")
    h5["openPMDextension"] = fstr("BeamPhysics;SpeciesType")
    h5["particlesPath"] = fstr(particlesPath)
end

"""
    pmd_field_init(h5, externalFieldPath="/ExternalFieldPath/%T/")

Root attribute initialization for an openPMD-beamphysics External Field Mesh

h5 should be the root of the file.
"""
function pmd_field_init(h5, externalFieldPath="/ExternalFieldPath/%T/")
    h5["dataType"] = fstr("openPMD")
    h5["openPMD"] = fstr("2.0.0")
    h5["openPMDextension"] = fstr("BeamPhysics")
    h5["externalFieldPath"] = fstr(externalFieldPath)
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
function write_pmd_bunch(h5, data; name=nothing)
    g = isnothing(name) ? h5 : create_group(h5, name)

    # Write into species group
    species = data["species"]
    g = create_group(g, species)

    # Attributes
    g["speciesType"] = fstr(species)
    g["numParticles"] = data["n_particle"]
    g["totalCharge"] = data["charge"]
    g["chargeUnitSI"] = 1.0

    # Required Datasets
    for key in ["x", "px", "y", "py", "z", "pz", "t", "status", "weight"]
        # Get full name, write data
        g2_name = COMPONENT_FROM_ALIAS[key]

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
    write_pmd_field(h5, data; name=nothing)

Data is a dict with:
    attrs: flat dict of attributes.
    components: flat dict of components

See inverse routine:
    .readers.load_field_data
"""
function write_pmd_field(h5, data; name=nothing)
    g = isnothing(name) ? h5 : create_group(h5, name)

    # Validate attrs
    attrs, other = load_field_attrs(data["attrs"])

    # Encode for writing
    attrs = encode_attrs(attrs)

    # Write attributes
    for (k, v) in attrs
        g[k] = v
    end

    # All other attributes (don't change them)
    for (k, v) in other
        g[k] = v
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
        g["value"] = dat0
        g["shape"] = size(data)
    else
        h5[name] = data
        g = h5[name]
    end

    if !isnothing(unit)
        g["unitSI"] = unit.unitSI
        g["unitDimension"] = unit.unitDimension
        g["unitSymbol"] = fstr(unit.unitSymbol)
    end

    return g
end