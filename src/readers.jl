"""
Utilities for reading particle and field data from HDF5 files in openPMD format.
"""

# -----------------------------------------
# Records, components, units

const PARTICLE_RECORD_COMPONENTS = Dict{String, Union{Nothing, Vector{String}}}(
    "branchIndex" => nothing,
    "chargeState" => nothing,
    "electricField" => ["x", "y", "z"],
    "elementIndex" => nothing,
    "magneticField" => ["x", "y", "z"],
    "locationInElement" => nothing,
    "momentum" => ["x", "y", "z"],
    "momentumOffset" => ["x", "y", "z"],
    "photonPolarizationAmplitude" => ["x", "y"],
    "photonPolarizationPhase" => ["x", "y"],
    "sPosition" => nothing,
    "totalMomentum" => nothing,
    "totalMomentumOffset" => nothing,
    "particleStatus" => nothing,
    "pathLength" => nothing,
    "position" => ["x", "y", "z"],
    "positionOffset" => ["x", "y", "z"],
    "spin" => ["x", "y", "z", "theta", "phi", "psi"],
    "time" => nothing,
    "timeOffset" => nothing,
    "velocity" => ["x", "y", "z"],
    "weight" => nothing
)

const FIELD_RECORD_COMPONENTS = Dict{String, Vector{String}}(
    "electricField" => ["x", "y", "z", "r", "theta"],
    "magneticField" => ["x", "y", "z", "r", "theta"]
)

# Expected unit dimensions for particle and field records
const EXPECTED_RECORD_UNIT_DIMENSION = Dict{String, NTuple{7, Int}}(
    "branchIndex" => dimension("1"),
    "chargeState" => dimension("1"),
    "electricField" => dimension("electric_field"),
    "magneticField" => dimension("magnetic_field"),
    "elementIndex" => dimension("1"),
    "locationInElement" => dimension("1"),
    "momentum" => dimension("momentum"),
    "momentumOffset" => dimension("momentum"),
    "photonPolarizationAmplitude" => dimension("electric_field"),
    "photonPolarizationPhase" => dimension("1"),
    "sPosition" => dimension("length"),
    "totalMomentum" => dimension("momentum"),
    "totalMomentumOffset" => dimension("momentum"),
    "particleStatus" => dimension("1"),
    "pathLength" => dimension("length"),
    "position" => dimension("length"),
    "positionOffset" => dimension("length"),
    "spin" => dimension("1"),
    "time" => dimension("time"),
    "timeOffset" => dimension("time"),
    "velocity" => dimension("velocity"),
    "weight" => dimension("charge")
)

# Convenient aliases for components
const COMPONENT_FROM_ALIAS = Dict{String, String}(
    "t" => "time",
    "weight" => "weight",
    "status" => "particleStatus"
)

# Add position, momentum, field aliases
for (g, prefix) in zip(["position", "momentum", "electricField", "magneticField"], ["", "p", "E", "B"])
    for c in ["x", "y", "z", "r", "theta"]
        alias = prefix * c
        COMPONENT_FROM_ALIAS[alias] = g * "/" * c
    end
end

# Inverse mapping
const COMPONENT_ALIAS = Dict(v => k for (k, v) in COMPONENT_FROM_ALIAS)

"""
    particle_paths(h5, key="particlesPath")

Uses the basePath and particlesPath to find where openPMD particles should be.
"""
function particle_paths(h5, key="particlesPath")
    basePath = String(h5.attrs["basePath"])
    particlesPath = String(h5.attrs[key])

    if !occursin("%T", basePath)
        return [basePath * particlesPath]
    end
    
    path1, path2 = split(basePath, "%T")
    tlist = keys(h5[path1])
    return [path1 * t * path2 * particlesPath for t in tlist]
end

"""
    field_paths(h5, key="externalFieldPath")

Looks for the External Fields.
"""
function field_paths(h5, key="externalFieldPath")
    if !haskey(h5.attrs, key)
        return String[]
    end

    fpath = String(h5.attrs[key])

    if !occursin("%T", fpath)
        return [fpath]
    end

    path1 = split(fpath, "%T")[1]
    tlist = keys(h5[path1])
    return [path1 * t for t in tlist]
end

"""
    is_constant_component(h5)

Constant record component should have 'value' and 'shape'.
"""
function is_constant_component(h5)
    return haskey(h5.attrs, "value") && haskey(h5.attrs, "shape")
end

"""
    constant_component_value(h5)

Get value from constant component.
"""
function constant_component_value(h5)
    unitSI = h5.attrs["unitSI"]
    val = h5.attrs["value"]
    return unitSI == 1.0 ? val : val * unitSI
end

"""
    component_unit_dimension(h5)

Return the unit dimension tuple.
"""
function component_unit_dimension(h5)
    return Tuple(h5.attrs["unitDimension"])
end

"""
    is_legacy_fortran_data_ordering(component_data_attrs)

Check if data is in legacy Fortran ordering.
"""
function is_legacy_fortran_data_ordering(component_data_attrs)
    if haskey(component_data_attrs, "gridDataOrder")
        @warn "Legacy gridDataOrder in component. Please remove and use axisLabels at the group level."
        if decode_attr(component_data_attrs["gridDataOrder"]) == "F"
            return true
        end
    end
    return false
end

"""
    component_data(h5, slice=:, unit_factor=1, axis_labels=nothing)

Returns an array from an h5 component.

Parameters:
- h5: HDF5.Dataset or HDF5.Group
- slice: Slice or tuple of slices to retrieve parts of the array
- unit_factor: Additional factor to convert from SI units to output units
- axis_labels: Required for multidimensional arrays. Supported options:
  * ("z", "y", "x")
  * ("z", "theta", "r")
  * ("x", "y", "z")
  * ("r", "theta", "z")
"""
function component_data(h5, slice=:, unit_factor=1, axis_labels=nothing)
    # Look for unitSI factor
    factor = haskey(h5.attrs, "unitSI") ? h5.attrs["unitSI"] : 1.0

    # Additional conversion factor
    if unit_factor != 1
        factor *= unit_factor
    end

    if is_constant_component(h5)
        dat = fill(h5.attrs["value"], h5.attrs["shape"])[slice]
    # Check multidimensional for data ordering, convert to 'x', 'y', 'z' order
    elseif ndims(h5) > 1
        if isnothing(axis_labels)
            throw(ArgumentError("axis_labels required for multidimensional arrays"))
        end

        # Reorder to x, y, z
        if axis_labels in [("z", "y", "x"), ("z", "theta", "r")]
            if slice isa Tuple
                # Need to transpose the slice ordering
                slice = reverse(slice)
            end

            # Retrieve dataset and transpose for C order
            dat = h5[slice]
            dat = permutedims(dat)
        elseif axis_labels in [("x", "y", "z"), ("r", "theta", "z")]
            dat = h5[slice]
        else
            # C-order
            dat = h5[slice]
        end
    else
        # 1-D array
        dat = h5[slice]
    end

    if factor != 1
        dat .*= factor
    end

    return dat
end

"""
    offset_component_name(component_name)

Many components can also have an offset, as in:
    position/x
    positionOffset/c
Return the appropriate name.
"""
function offset_component_name(component_name)
    x = split(component_name, "/")
    if length(x) == 1
        return x[1] * "Offset"
    else
        return x[1] * "Offset/" * x[2]
    end
end

"""
    particle_array(h5, component, slice=:, include_offset=true)

Main routine to return particle arrays in fixed units.
All units are SI except momentum, which will be in eV/c.
"""
function particle_array(h5, component, slice=:, include_offset=true)
    # Handle aliases
    if haskey(COMPONENT_FROM_ALIAS, component)
        component = COMPONENT_FROM_ALIAS[component]
    end

    if component in ["momentum/x", "momentum/y", "momentum/z"]
        unit_factor = c_light / e_charge  # convert J/(m/s) to eV/c
    else
        unit_factor = 1.0
    end

    # Get data
    dat = component_data(h5[component], slice=slice, unit_factor=unit_factor)

    # Look for offset component
    ocomponent = offset_component_name(component)
    if include_offset && haskey(h5, ocomponent)
        offset = component_data(h5[ocomponent], slice=slice, unit_factor=unit_factor)
        dat .+= offset
    end

    return dat
end

"""
    all_components(h5)

Look for possible components in a particle group.
"""
function all_components(h5)
    components = String[]
    for record_name in keys(h5)
        if !haskey(PARTICLE_RECORD_COMPONENTS, record_name)
            continue
        end

        # Look for components
        possible_components = PARTICLE_RECORD_COMPONENTS[record_name]

        if isnothing(possible_components)
            # Record is a component
            push!(components, record_name)
        else
            g = h5[record_name]
            for cname in possible_components
                if haskey(g, cname)
                    push!(components, record_name * "/" * cname)
                end
            end
        end
    end
    return components
end

"""
    component_str(particle_group, name)

Informational string from a component in a particle group (h5).
"""
function component_str(particle_group, name)
    g = particle_group[name]
    record_name = split(name, "/")[1]
    expected_dimension = EXPECTED_RECORD_UNIT_DIMENSION[record_name]
    this_dimension = component_unit_dimension(g)
    dname = dimension_name(this_dimension)
    symbol = SI_SYMBOL[dname]

    s = name * " "

    if is_constant_component(g)
        val = constant_component_value(g)
        shape = g.attrs["shape"]
        s *= "[constant $val with shape $shape]"
    else
        s *= "[" * string(length(g)) * " items]"
    end

    if symbol != "1"
        s *= " is a $dname with units: $symbol"
    end

    if expected_dimension != this_dimension
        s *= ", but expected units: " * SI_SYMBOL[dimension_name(expected_dimension)]
    end

    return s
end

# ----------------------------------
# Fields

const REQUIRED_FIELD_ATTRS = [
    # strings
    "eleAnchorPt",
    "gridGeometry",
    "axisLabels",
    # reals and ints
    "gridLowerBound",
    "gridOriginOffset",
    "gridSpacing",
    "gridSize",
    "harmonic"
]

# Dict with options
const OPTIONAL_FIELD_ATTRS = Dict{String, Any}(
    "name" => nothing,
    "gridCurvatureRadius" => nothing,
    "fundamentalFrequency" => 0,
    "RFphase" => 0,
    "fieldScale" => 1.0,
    "masterParameter" => nothing
)

"""
    load_field_attrs(attr, verbose=false)

Loads FieldMesh required and optional attributes from a dict_like object.
Non-standard attributes will be collected in an 'other' dict.

Returns tuple: (attrs, other)
"""
function load_field_attrs(attr, verbose=false)
    # Get all attrs. Will pop.
    a = Dict(attr)

    attrs = Dict{String, Any}()
    other = Dict{String, Any}()

    # Required
    for k in REQUIRED_FIELD_ATTRS
        attrs[k] = pop!(a, k)
    end

    # Optional, filling in some defaults
    for (k, v) in OPTIONAL_FIELD_ATTRS
        if haskey(a, k)
            attrs[k] = pop!(a, k)
        elseif !isnothing(v)
            attrs[k] = v
        end
    end

    # Collect other
    for (k, v) in a
        other[k] = v
        if verbose
            println("Nonstandard attr: ", k, " ", v)
        end
    end

    # Decode
    attrs = decode_attrs(attrs)

    return attrs, other
end