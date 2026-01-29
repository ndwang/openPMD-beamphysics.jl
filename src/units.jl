"""
Simple units functionality for the openPMD beamphysics records.

This module provides utilities for handling units in openPMD beamphysics records, including:
- Conversion between Unitful.jl units and openPMD 7-tuple dimension format
- Unit conversion and scaling
- HDF5 I/O with unit support
- SI prefix handling

For more advanced units functionality, use Unitful.jl directly.
"""

# Map Unitful dimensions to openPMD 7-tuple format
const DIMENSION_MAP = Dict{Type, Int}(
    Unitful.Dimension{:Length} => 1,
    Unitful.Dimension{:Mass} => 2,
    Unitful.Dimension{:Time} => 3,
    Unitful.Dimension{:Current} => 4,
    Unitful.Dimension{:Temperature} => 5,
    Unitful.Dimension{:Amount} => 6,
    Unitful.Dimension{:Luminosity} => 7
)

# Base dimensions for openPMD
const DIMENSION = Dict{String, NTuple{7, Int}}(
    "1" => (0, 0, 0, 0, 0, 0, 0),
    # Base units
    "length" => (1, 0, 0, 0, 0, 0, 0),
    "mass" => (0, 1, 0, 0, 0, 0, 0),
    "time" => (0, 0, 1, 0, 0, 0, 0),
    "current" => (0, 0, 0, 1, 0, 0, 0),
    "temperature" => (0, 0, 0, 0, 1, 0, 0),
    "mol" => (0, 0, 0, 0, 0, 1, 0),
    "luminous" => (0, 0, 0, 0, 0, 0, 1),
    #
    "charge" => (0, 0, 1, 1, 0, 0, 0),
    "electric_field" => (1, 1, -3, -1, 0, 0, 0),
    "electric_potential" => (1, 2, -3, -1, 0, 0, 0),
    "magnetic_field" => (0, 1, -2, -1, 0, 0, 0),
    "velocity" => (1, 0, -1, 0, 0, 0, 0),
    "energy" => (2, 1, -2, 0, 0, 0, 0),
    "momentum" => (1, 1, -1, 0, 0, 0, 0)
)

# Inverse mapping
const DIMENSION_NAME = Dict{NTuple{7, Int}, String}(v => k for (k, v) in DIMENSION)

# SI symbols
const SI_SYMBOL = Dict{String, String}(
    "1" => "1",
    "length" => "m",
    "mass" => "kg",
    "time" => "s",
    "current" => "A",
    "temperature" => "K",
    "mol" => "mol",
    "luminous" => "cd",
    "charge" => "C",
    "electric_field" => "V/m",
    "electric_potential" => "V",
    "velocity" => "m/s",
    "energy" => "J",
    "momentum" => "kg⋅m/s",
    "magnetic_field" => "T"
)

# Inverse mapping
const SI_NAME = Dict{String, String}(v => k for (k, v) in SI_SYMBOL)

# Known units mapping to Unitful units
const KNOWN_UNITS = Dict{String, Unitful.FreeUnits}(
    "1" => NoUnits,
    "degree" => u"°",
    "rad" => u"rad",
    "m" => u"m",
    "kg" => u"kg",
    "g" => u"g",
    "s" => u"s",
    "A" => u"A",
    "K" => u"K",
    "mol" => u"mol",
    "cd" => u"cd",
    "C" => u"C",
    "V/m" => u"V/m",
    "V" => u"V",
    "c_light" => u"c",
    "m/s" => u"m/s",
    "eV" => u"eV",
    "J" => u"J",
    "eV/c" => u"eV/c",
    "W" => u"W",
    "W/rad^2" => u"W/rad^2",
    "W/m^2" => u"W/m^2",
    "T" => u"T"
)

# Prefix dictionaries
const PREFIX_FACTOR = Dict{String, Float64}(
    "yocto-" => 1e-24,
    "zepto-" => 1e-21,
    "atto-" => 1e-18,
    "femto-" => 1e-15,
    "pico-" => 1e-12,
    "nano-" => 1e-9,
    "micro-" => 1e-6,
    "milli-" => 1e-3,
    "centi-" => 1e-2,
    "deci-" => 1e-1,
    "deca-" => 1e1,
    "hecto-" => 1e2,
    "kilo-" => 1e3,
    "mega-" => 1e6,
    "giga-" => 1e9,
    "tera-" => 1e12,
    "peta-" => 1e15,
    "exa-" => 1e18,
    "zetta-" => 1e21,
    "yotta-" => 1e24
)

# Inverse
const PREFIX = Dict{Float64, String}(v => k for (k, v) in PREFIX_FACTOR)

const SHORT_PREFIX_FACTOR = Dict{String, Float64}(
    "y" => 1e-24,
    "z" => 1e-21,
    "a" => 1e-18,
    "f" => 1e-15,
    "p" => 1e-12,
    "n" => 1e-9,
    "µ" => 1e-6,
    "m" => 1e-3,
    "c" => 1e-2,
    "d" => 1e-1,
    "" => 1,
    "da" => 1e1,
    "h" => 1e2,
    "k" => 1e3,
    "M" => 1e6,
    "G" => 1e9,
    "T" => 1e12,
    "P" => 1e15,
    "E" => 1e18,
    "Z" => 1e21,
    "Y" => 1e24
)

# Inverse
const SHORT_PREFIX = Dict{Float64, String}(v => k for (k, v) in SHORT_PREFIX_FACTOR)

"""
    dimension(name::String) -> NTuple{7, Int}

Get the dimension tuple for a named quantity.

The dimension tuple follows the openPMD 7-tuple format:
1. Length
2. Mass
3. Time
4. Current
5. Temperature
6. Amount
7. Luminosity
"""
function dimension(name::String)
    haskey(DIMENSION, name) || error("Invalid unit dimension string: $name. Valid options are: $(join(keys(DIMENSION), ", "))")
    return DIMENSION[name]
end

"""
    UnitfulToOpenPMD(u::Unitful.FreeUnits) -> NTuple{7, Int}

Convert a Unitful.jl unit to openPMD 7-tuple dimension format.
"""
function UnitfulToOpenPMD(u::Unitful.FreeUnits)
    # Get the dimension from Unitful
    # dim is a tuple of dimension objects
    dim = typeof(Unitful.dimension(u)).parameters[1]
    # Extract the power of each dimension
    dim_array = zeros(Int, 7)
    for d in dim
        idx = get(DIMENSION_MAP, typeof(d)) do
            error("Unknown dimension type: $(typeof(d))")
        end
        dim_array[idx] = d.power
    end
    return Tuple(dim_array)
end

"""
    dimension_name(dim_array::NTuple{7, Int}) -> String

Get the name of a dimension tuple.
"""
function dimension_name(dim_array::NTuple{7, Int})
    haskey(DIMENSION_NAME, dim_array) || error("Unknown dimension tuple: $dim_array")
    return DIMENSION_NAME[dim_array]
end

"""
    dimension_name(u::Unitful.FreeUnits) -> String

Get the name of a Unitful unit's dimensions.
"""
function dimension_name(u::Unitful.FreeUnits)
    return dimension_name(UnitfulToOpenPMD(u))
end

"""
    SI_conversion_factor(unit::Unitful.Units) -> Float64

Returns the conversion factor to convert from the given unit to SI units.
"""
function SI_conversion_factor(unit::Unitful.Units)
    # 1. Create a quantity with magnitude 1 of the given unit
    quantity_one = 1.0 * unit
    # 2. Convert to preferred SI units (Unitful's default preferred system is SI)
    si_quantity = upreferred(quantity_one)
    # 3. Strip the units to get the numerical factor
    factor = ustrip(si_quantity)

    return factor
end

"""
    unit_conversion_factor(src_unit::Unitful.FreeUnits, dst_unit::Unitful.FreeUnits) -> Float64

Calculate the conversion factor between two units.
"""
function unit_conversion_factor(src_unit::Unitful.FreeUnits, dst_unit::Unitful.FreeUnits)
    # 1. Create a quantity with magnitude 1 of src_unit
    quantity_one = 1.0 * src_unit
    # 2. Convert to dst_unit
    si_quantity = uconvert(dst_unit, quantity_one)
    # 3. Strip the units to get the numerical factor
    factor = ustrip(si_quantity)

    return factor
end

"""
    nice_scale_prefix(scale::Real) -> Tuple{Float64, String}

Returns a nice factor and an SI prefix string.

# Arguments
- `scale::Real`: The value to scale

# Returns
- scaling::Float64: The nice scaling factor
- prefix::String: The SI prefix string corresponding to the scale
"""
function nice_scale_prefix(scale::Real)
    if scale == 0
        return 1.0, ""
    end

    p10 = log10(abs(scale))

    if p10 < -24  # Limits of SI prefixes
        f = 1e-24
    elseif p10 > 24
        f = 1e24
    elseif p10 < -1.5 || p10 > 2
        f = 10.0^(floor(Int, p10/3) * 3)
    else
        f = 1.0
    end

    return f, SHORT_PREFIX[f]
end

"""
    nice_array(a::AbstractArray) -> Tuple{AbstractArray, Float64, String}

Scale an input array and return the scaled array, the scaling factor, and the
corresponding unit prefix.

# Arguments
- `a::AbstractArray`: The input array to scale

# Returns
- scaled_array::AbstractArray: The scaled array
- scaling::Float64: The scaling factor
- prefix::String: The SI prefix string corresponding to the scale
"""
function nice_array(a::AbstractArray)
    if length(a) == 1
        x = a[1]
    else
        x = maximum((maximum(a)-minimum(a), abs(mean(a))))  # Account for tiny spread
    end

    fac, prefix = nice_scale_prefix(x)
    return a ./ fac, fac, prefix
end

"""
    nice_array(a) -> Tuple{typeof(a), Float64, String}

Scale an input scalar and return the scaled scalar, the scaling factor, and the
corresponding unit prefix.

# Arguments
- `a`: The input scalar to scale

# Returns
- scaled_scalar: The scaled scalar
- scaling::Float64: The scaling factor
- prefix::String: The SI prefix string corresponding to the scale
"""
function nice_array(a)
    fac, prefix = nice_scale_prefix(a)
    return a / fac, fac, prefix
end


"""
    limits(x::AbstractArray) -> Tuple

Get the limits of an array.
"""
function limits(x::AbstractArray)
    return minimum(x), maximum(x)
end

"""
    plottable_array(x::AbstractArray; lim::Union{Nothing, Tuple{Real, Real}}=nothing) -> Tuple{AbstractArray, Float64, String, Real, Real}

Scale an input array and return the scaled array, the scaling factor, the
corresponding unit prefix, and the limits of the array.

# Arguments
- `x::AbstractArray`: The input array to scale
- `lim::Union{Nothing, Tuple{Real, Real}}`: The limits of the array

# Returns
- scaled_array::AbstractArray: The scaled array
- scaling::Float64: The scaling factor
- prefix::String: The SI prefix string corresponding to the scale
- xmin::Real: The minimum value of the array
- xmax::Real: The maximum value of the array
"""
function plottable_array(x::AbstractArray; lim::Union{Nothing, Tuple{Real, Real}}=nothing)
    xmin, xmax = lim === nothing ? limits(x) : lim

    _, factor, p1 = nice_array([xmin, xmax])

    return x ./ factor, factor, p1, xmin, xmax
end

"""
    write_unit_h5(h5, u::Unitful.FreeUnits)

Writes a Unitful unit to an HDF5 handle in openPMD format.

# Arguments
- `h5`: HDF5 handle to write to
- `u::Unitful.FreeUnits`: The unit to write
"""
function write_unit_h5(h5, u::Unitful.FreeUnits)
    h5["unitSI"] = SI_conversion_factor(u)
    h5["unitDimension"] = collect(dimension(dimension_name(u)))  # HDF5.jl only accepts arrays for attributes
    h5["unitSymbol"] = string(u) # TODO: uparse(string(u)) does not work on eV/c
end

"""
    read_unit_h5(h5) -> Unitful.FreeUnits

Reads unit data from an HDF5 handle and returns a Unitful FreeUnits object.

# Arguments
- `h5`: HDF5 handle to read from

# Returns
- `Unitful.FreeUnits`: The unit read from the HDF5 handle
"""
function read_unit_h5(h5)
    unitSymbol = read(h5["unitSymbol"])
    return Unitful.uparse(unitSymbol)
end

"""
    read_dataset_and_unit_h5(h5)

Reads a dataset that has openPMD unit attributes.

# Arguments
- `h5`: HDF5 handle to read from

# Returns
- data: The data array
- unit::Unitful.FreeUnits: The unit of the data
"""
function read_dataset_and_unit_h5(h5)
    u = read_unit_h5(h5)
    return Array(h5), u
end

"""
    read_dataset_and_unit_h5(h5, expected_unit::String; convert::Bool=true)

Reads a dataset that has openPMD unit attributes and converts it to the expected unit.

# Arguments
- `h5`: HDF5 handle to read from
- `expected_unit::String`: The expected unit
- `convert::Bool`: Whether to convert the data to the expected unit (default: true)

# Returns
- data: The data array
- unit::Unitful.FreeUnits: The unit of the data

# Throws
- `ErrorException`: If the units are incompatible when expected_unit is provided
"""
function read_dataset_and_unit_h5(h5, expected_unit::String; convert::Bool=true)
    # Read the unit that is there
    u = read_unit_h5(h5)
    expected_unit = Unitful.uparse(expected_unit)

    # Check dimensions
    if Unitful.dimension(u) != Unitful.dimension(expected_unit)
        error("Unit dimensions do not match: $(u) vs $(expected_unit)")
    end

    if convert
        fac = unit_conversion_factor(u, expected_unit)
        return fac * Array(h5), expected_unit
    else
        return Array(h5), u
    end
end

"""
    read_dataset_and_unit_h5(h5, expected_unit::Unitful.FreeUnits; convert::Bool=true)

Reads a dataset that has openPMD unit attributes and converts it to the expected unit.

# Arguments
- `h5`: HDF5 handle to read from
- `expected_unit::Unitful.FreeUnits`: The expected unit
- `convert::Bool`: Whether to convert the data to the expected unit (default: true)

# Returns
- data: The data array
- unit::Unitful.FreeUnits: The unit of the data

# Throws
- `ErrorException`: If the units are incompatible when expected_unit is provided
"""
function read_dataset_and_unit_h5(h5, expected_unit::Unitful.FreeUnits; convert::Bool=true)
    # Read the unit that is there
    u = read_unit_h5(h5)

    # Check dimensions
    if Unitful.dimension(u) != Unitful.dimension(expected_unit)
        error("Unit dimensions do not match: $(u) vs $(expected_unit)")
    end

    if convert
        fac = unit_conversion_factor(u, expected_unit)
        return fac * Array(h5), expected_unit
    else
        return Array(h5), u
    end
end

"""
    write_dataset_and_unit_h5(h5, name::String, data; unit=nothing)

Writes data and unit to an HDF5 dataset.

# Arguments
- `h5`: HDF5 handle to write to
- `name::String`: Name of the dataset
- `data`: The data to write
- `unit`: The unit of the data (default: nothing)
"""
function write_dataset_and_unit_h5(h5, name::String, data; unit=nothing)
    h5[name] = data
    if unit !== nothing
        write_unit_h5(h5[name], unit)
    end
end

# ParticleGroup units
const PARTICLEGROUP_UNITS = Dict{String, Unitful.FreeUnits}(
    # Basic quantities
    "n_particle" => NoUnits,
    "status" => NoUnits,
    "id" => NoUnits,
    "n_alive" => NoUnits,
    "n_dead" => NoUnits,
    
    # Time
    "t" => u"s",
    "z/c" => u"s",
    
    # Energy and mass
    "energy" => u"eV",
    "kinetic_energy" => u"eV",
    "mass" => u"eV",
    "higher_order_energy_spread" => u"eV",
    "higher_order_energy" => u"eV",
    
    # Momentum
    "px" => u"eV/c",
    "py" => u"eV/c",
    "pz" => u"eV/c",
    "p" => u"eV/c",
    "pr" => u"eV/c",
    "ptheta" => u"eV/c",
    
    # Position
    "x" => u"m",
    "y" => u"m",
    "z" => u"m",
    "r" => u"m",
    "Jx" => u"m",
    "Jy" => u"m",
    
    # Dimensionless quantities
    "beta" => NoUnits,
    "beta_x" => NoUnits,
    "beta_y" => NoUnits,
    "beta_z" => NoUnits,
    "gamma" => NoUnits,
    "bunching" => NoUnits,
    
    # Angles
    "theta" => u"rad",
    "bunching_phase" => u"rad",
    
    # Charge
    "charge" => u"C",
    "species_charge" => u"C",
    "weight" => u"C",
    "average_current" => u"A",
    
    # Emittance
    "norm_emit_x" => u"m",
    "norm_emit_y" => u"m",
    "norm_emit_4d" => u"m^2",
    
    # Angular momentum
    "Lz" => u"m*eV/c",
    
    # Angles
    "xp" => u"rad",
    "yp" => u"rad",

    # Normalized coordinates
    "x_bar" => u"m^(1/2)",
    "px_bar" => u"m^(1/2)",
    "y_bar" => u"m^(1/2)",
    "py_bar" => u"m^(1/2)"
)

# Add field components
for component in ["", "x", "y", "z", "theta", "r"]
    PARTICLEGROUP_UNITS["E$(component)"] = u"V/m"
    PARTICLEGROUP_UNITS["B$(component)"] = u"T"
end

# Add Twiss parameters
for plane in ("x", "y")
    PARTICLEGROUP_UNITS["twiss_alpha_$(plane)"] = NoUnits
    PARTICLEGROUP_UNITS["twiss_etap_$(plane)"] = NoUnits
    PARTICLEGROUP_UNITS["twiss_beta_$(plane)"] = u"m"
    PARTICLEGROUP_UNITS["twiss_eta_$(plane)"] = u"m"
    PARTICLEGROUP_UNITS["twiss_emit_$(plane)"] = u"m"
    PARTICLEGROUP_UNITS["twiss_norm_emit_$(plane)"] = u"m"
    PARTICLEGROUP_UNITS["twiss_gamma_$(plane)"] = u"m^-1"
end

"""
    pg_units(key::String) -> Unitful.FreeUnits

Returns the unit for any attribute in a ParticleGroup.

# Arguments
- `key::String`: The attribute key to get the unit for

# Returns
- `Unitful.FreeUnits`: The unit for the attribute

# Throws
- `ErrorException`: If no unit is known for the key
"""
function pg_units(key::String)
    # Basic cases
    if haskey(PARTICLEGROUP_UNITS, key)
        return PARTICLEGROUP_UNITS[key]
    end

    # Operators
    for prefix in ["sigma_", "mean_", "min_", "max_", "ptp_", "delta_"]
        if startswith(key, prefix)
            nkey = key[length(prefix)+1:end]
            return pg_units(nkey)
        end
    end

    # Covariance
    if startswith(key, "cov_")
        subkeys = split(key[5:end], "__")
        length(subkeys) == 2 || error("Invalid covariance key format: $key")
        unit0 = PARTICLEGROUP_UNITS[subkeys[1]]
        unit1 = PARTICLEGROUP_UNITS[subkeys[2]]
        return unit0 * unit1
    end

    # Fields
    if startswith(key, "electricField")
        return u"V/m"
    end
    if startswith(key, "magneticField")
        return u"T"
    end
    if startswith(key, "bunching_phase")
        return u"rad"
    end
    if startswith(key, "bunching")
        return NoUnits
    end

    error("No known unit for: $key")
end