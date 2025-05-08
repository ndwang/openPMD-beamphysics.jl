"""
Simple units functionality for the openPMD beamphysics records.

For more advanced units, use a package like Unitful.jl:
    https://github.com/PainterQubits/Unitful.jl
"""

# Constants
# TODO: PhysicalConstants.jl uses Unitful.Quantity
#       Need to modify this package to work with Quantity
#       or convert to Float64

const c_light = Float64(2.99792458e8)
const e_charge = Float64(1.602176634e-19)
const m_e = Float64(9.1093837015e-31)
const m_p = Float64(1.67262192369e-27)
const m_mu = Float64(1.83615267e-28)

const mec2 = m_e * c_light ^2 / e_charge
const mpc2 = m_p * c_light ^2 / e_charge
const mmc2 = m_mu * c_light ^2 / e_charge
const mhmc2 = mpc2 + mec2 * 2
const mH2pc2 = 2 * mpc2 + mec2

const mu_0 = Float64(1.25663706212e-6)

# Type aliases
const Limit = Tuple{Union{Float64, Nothing}, Union{Float64, Nothing}}
const Dimension = NTuple{7, Int}

"""
    PMDUnit

OpenPMD representation of a unit.

# Fields
- `unitSymbol::String`: Native units name
- `unitSI::Float64`: Conversion factor to the corresponding SI unit
- `unitDimension::Dimension`: 7-tuple of SI Base Exponents

# Examples
```julia
# Define that an eV is 1.602176634e-19 of base units m^2 kg/s^2, which is a Joule (J)
PMDUnit("eV", 1.602176634e-19, (2, 1, -2, 0, 0, 0, 0))

# If unitSI=0 (default), PMDUnit may be initialized with a known symbol
PMDUnit("T")
```
"""
struct PMDUnit
    unitSymbol::String
    unitSI::Float64
    unitDimension::Dimension

    function PMDUnit(unitSymbol::String="", unitSI::Real=0, unitDimension::Union{String, Dimension}=(0, 0, 0, 0, 0, 0, 0))
        # Allow to return an internally known unit
        if unitSI == 0
            if !haskey(KNOWN_UNIT, unitSymbol)
                throw(ArgumentError("Unknown unitSymbol: $unitSymbol"))
            end

            # Copy internals
            u = KNOWN_UNIT[unitSymbol]
            unitSI = u.unitSI
            unitDimension = u.unitDimension
        end

        if unitDimension isa String
            unitDimension = dimension(unitDimension)
        else
            unitDimension = make_dimension(unitDimension)
        end

        new(unitSymbol, Float64(unitSI), unitDimension)
    end
end

# Base methods for PMDUnit
Base.:(==)(u1::PMDUnit, u2::PMDUnit) = u1.unitSymbol == u2.unitSymbol && 
                                      u1.unitSI == u2.unitSI && 
                                      u1.unitDimension == u2.unitDimension

Base.:(*)(u1::PMDUnit, u2::PMDUnit) = multiply_units(u1, u2)
Base.:(/)(u1::PMDUnit, u2::PMDUnit) = divide_units(u1, u2)

Base.show(io::IO, u::PMDUnit) = print(io, u.unitSymbol)
Base.show(io::IO, ::MIME"text/plain", u::PMDUnit) = 
    print(io, "PMDUnit('$(u.unitSymbol)', $(u.unitSI), $(u.unitDimension))")

"""
    is_dimensionless(u::PMDUnit)

Checks if the unit is dimensionless
"""
is_dimensionless(u::PMDUnit) = u.unitDimension == (0, 0, 0, 0, 0, 0, 0)

"""
    is_identity(u::PMDUnit)

Checks if the unit is equivalent to 1
"""
is_identity(u::PMDUnit) = u.unitSI == 1 && u.unitDimension == (0, 0, 0, 0, 0, 0, 0)

"""
    multiply_units(u1::PMDUnit, u2::PMDUnit)

Multiplies two PMDUnit symbols.
"""
function multiply_units(u1::PMDUnit, u2::PMDUnit)
    if is_identity(u1)
        return u2
    end
    if is_identity(u2)
        return u1
    end

    s1 = u1.unitSymbol
    s2 = u2.unitSymbol
    symbol = s1 == s2 ? "($s1)^2" : s1 * "*" * s2
    dim = Tuple(sum(x) for x in zip(u1.unitDimension, u2.unitDimension))
    unitSI = u1.unitSI * u2.unitSI

    return PMDUnit(unitSymbol=symbol, unitSI=unitSI, unitDimension=dim)
end

"""
    divide_units(u1::PMDUnit, u2::PMDUnit)

Divides two PMDUnit symbols: u1 / u2
"""
function divide_units(u1::PMDUnit, u2::PMDUnit)
    if is_identity(u2)
        return u1
    end

    s1 = u1.unitSymbol
    s2 = u2.unitSymbol
    symbol = s1 == s2 ? "1" : s1 * "/" * s2
    dim = make_dimension(a - b for (a, b) in zip(u1.unitDimension, u2.unitDimension))
    unitSI = u1.unitSI / u2.unitSI

    return PMDUnit(unitSymbol=symbol, unitSI=unitSI, unitDimension=dim)
end

"""
    sqrt_unit(u::PMDUnit)

Returns the square root of a unit.
"""
function sqrt_unit(u::PMDUnit)
    symbol = u.unitSymbol
    if symbol ∉ ["", "1"]
        symbol = "\\sqrt{ $symbol }"
    end

    unitSI = sqrt(u.unitSI)
    dim = Tuple(x ÷ 2 for x in u.unitDimension)

    return PMDUnit(unitSymbol=symbol, unitSI=unitSI, unitDimension=dim)
end

# Dimension definitions
const DIMENSION = Dict{String, Dimension}(
    "1" => (0, 0, 0, 0, 0, 0, 0),
    # Base units
    "length" => (1, 0, 0, 0, 0, 0, 0),
    "mass" => (0, 1, 0, 0, 0, 0, 0),
    "time" => (0, 0, 1, 0, 0, 0, 0),
    "current" => (0, 0, 0, 1, 0, 0, 0),
    "temperture" => (0, 0, 0, 0, 1, 0, 0),
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

# Inverse
const DIMENSION_NAME = Dict{Dimension, String}(v => k for (k, v) in DIMENSION)

"""
    make_dimension(dim)

Convert a sequence of integers to a Dimension tuple.
"""
function make_dimension(dim)
    dim = Tuple(Int(d) for d in dim)
    if length(dim) != 7
        throw(ArgumentError("Dimensions must be 7 elements."))
    end
    return dim
end

"""
    dimension(name::String)

Get the dimension tuple for a given dimension name.
"""
function dimension(name::String)
    if !haskey(DIMENSION, name)
        options = join(keys(DIMENSION), ", ")
        throw(ArgumentError("Invalid unit dimension string: $name. Valid options are: $options"))
    end
    return DIMENSION[name]
end

"""
    dimension_name(dim_array::Dimension)

Get the dimension name for a given dimension tuple.
"""
dimension_name(dim_array::Dimension) = DIMENSION_NAME[make_dimension(dim_array)]

# SI symbols
const SI_SYMBOL = Dict{String, String}(
    "1" => "1",
    "length" => "m",
    "mass" => "kg",
    "time" => "s",
    "current" => "A",
    "temperture" => "K",
    "mol" => "mol",
    "luminous" => "cd",
    "charge" => "C",
    "electric_field" => "V/m",
    "electric_potential" => "V",
    "velocity" => "m/s",
    "energy" => "J",
    "momentum" => "kg*m/s",
    "magnetic_field" => "T"
)

# Inverse
const SI_NAME = Dict{String, String}(v => k for (k, v) in SI_SYMBOL)

# Known units
const KNOWN_UNIT = Dict{String, PMDUnit}(
    "1" => PMDUnit("", 1, "1"),
    "degree" => PMDUnit("degree", π/180, "1"),
    "rad" => PMDUnit("rad", 1, "1"),
    "m" => PMDUnit("m", 1, "length"),
    "kg" => PMDUnit("kg", 1, "mass"),
    "g" => PMDUnit("g", 0.001, "mass"),
    "s" => PMDUnit("s", 1, "time"),
    "A" => PMDUnit("A", 1, "current"),
    "K" => PMDUnit("K", 1, "temperture"),
    "mol" => PMDUnit("mol", 1, "mol"),
    "cd" => PMDUnit("cd", 1, "luminous"),
    "C" => PMDUnit("C", 1, "charge"),
    "charge_num" => PMDUnit("charge #", 1, "charge"),
    "V/m" => PMDUnit("V/m", 1, "electric_field"),
    "V" => PMDUnit("V", 1, "electric_potential"),
    "c_light" => PMDUnit("vel/c", c_light, "velocity"),
    "m/s" => PMDUnit("m/s", 1, "velocity"),
    "eV" => PMDUnit("eV", e_charge, "energy"),
    "J" => PMDUnit("J", 1, "energy"),
    "eV/c" => PMDUnit("eV/c", e_charge/c_light, "momentum"),
    "eV/m" => PMDUnit("eV/m", e_charge, (1, 1, -2, 0, 0, 0, 0)),
    "W" => PMDUnit("W", 1, (2, 1, -3, 0, 0, 0, 0)),
    "W/rad^2" => PMDUnit("W/rad^2", 1, (2, 1, -3, 0, 0, 0, 0)),
    "W/m^2" => PMDUnit("W/m^2", 1, (0, 1, -3, 0, 0, 0, 0)),
    "T" => PMDUnit("T", 1, "magnetic_field")
)

"""
    unit(symbol::String)

Returns a PMDUnit from a known symbol.

* and / are allowed between two known symbols.
"""
function unit(symbol::String)
    if haskey(KNOWN_UNIT, symbol)
        return KNOWN_UNIT[symbol]
    end

    if '*' in symbol || '/' in symbol
        parts = split(symbol, r"([*/])")
        result = KNOWN_UNIT[parts[1]]
        for (op, part) in zip(parts[2:2:end], parts[3:2:end])
            unit = KNOWN_UNIT[part]
            if op == "*"
                result = result * unit
            elseif op == "/"
                result = result / unit
            end
        end
        return result
    end

    throw(ArgumentError("Unknown unit symbol: $symbol"))
end

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
    nice_scale_prefix(scale::Real)

Returns a nice factor and an SI prefix string.
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
    nice_array(a::AbstractArray)

Scale an input array and return the scaled array, the scaling factor, and the
corresponding unit prefix.
"""
function nice_array(a::AbstractArray)
    if isscalar(a)
        x = a
    elseif length(a) == 1
        x = a[1]
    else
        x = max(peak2peak(a), abs(mean(a)))  # Account for tiny spread
    end

    fac, prefix = nice_scale_prefix(x)
    return a ./ fac, fac, prefix
end

"""
    plottable_array(x::AbstractArray; nice::Bool=true, lim::Union{Limit, Nothing}=nothing)

Similar to nice_array, but also considers limits for plotting.
"""
function plottable_array(x::AbstractArray; nice::Bool=true, lim::Union{Limit, Nothing}=nothing)
    x = collect(x)
    if !isnothing(lim)
        xmin = isnothing(lim[1]) ? minimum(x) : lim[1]
        xmax = isnothing(lim[2]) ? maximum(x) : lim[2]
    else
        xmin = minimum(x)
        xmax = maximum(x)
    end

    if nice
        _, factor, p1 = nice_array([xmin, xmax])
    else
        factor, p1 = 1.0, ""
    end

    return x ./ factor, factor, p1, xmin, xmax
end

# ParticleGroup units
const PARTICLEGROUP_UNITS = Dict{String, PMDUnit}()

# Initialize ParticleGroup units
for k in ["n_particle", "status", "id", "n_alive", "n_dead"]
    PARTICLEGROUP_UNITS[k] = unit("1")
end

for k in ["t", "z/c"]
    PARTICLEGROUP_UNITS[k] = unit("s")
end

for k in ["energy", "kinetic_energy", "mass", "higher_order_energy_spread", "higher_order_energy"]
    PARTICLEGROUP_UNITS[k] = unit("eV")
end

for k in ["px", "py", "pz", "p", "pr", "ptheta"]
    PARTICLEGROUP_UNITS[k] = unit("eV/c")
end

for k in ["x", "y", "z", "r", "Jx", "Jy"]
    PARTICLEGROUP_UNITS[k] = unit("m")
end

for k in ["beta", "beta_x", "beta_y", "beta_z", "gamma", "bunching"]
    PARTICLEGROUP_UNITS[k] = unit("1")
end

for k in ["theta", "bunching_phase"]
    PARTICLEGROUP_UNITS[k] = unit("rad")
end

for k in ["charge", "species_charge", "weight"]
    PARTICLEGROUP_UNITS[k] = unit("C")
end

for k in ["average_current"]
    PARTICLEGROUP_UNITS[k] = unit("A")
end

for k in ["norm_emit_x", "norm_emit_y"]
    PARTICLEGROUP_UNITS[k] = unit("m")
end

PARTICLEGROUP_UNITS["norm_emit_4d"] = multiply_units(unit("m"), unit("m"))
PARTICLEGROUP_UNITS["Lz"] = multiply_units(unit("m"), unit("eV/c"))

for k in ["xp", "yp"]
    PARTICLEGROUP_UNITS[k] = unit("rad")
end

for k in ["x_bar", "px_bar", "y_bar", "py_bar"]
    PARTICLEGROUP_UNITS[k] = sqrt_unit(unit("m"))
end

for component in ["", "x", "y", "z", "theta", "r"]
    PARTICLEGROUP_UNITS["E$component"] = unit("V/m")
    PARTICLEGROUP_UNITS["B$component"] = unit("T")
end

# Twiss
for plane in ("x", "y")
    for k in ("alpha", "etap")
        PARTICLEGROUP_UNITS["twiss_$(k)_$(plane)"] = unit("1")
    end
    for k in ("beta", "eta", "emit", "norm_emit")
        PARTICLEGROUP_UNITS["twiss_$(k)_$(plane)"] = unit("m")
    end
    for k in ("gamma",)
        PARTICLEGROUP_UNITS["twiss_$(k)_$(plane)"] = divide_units(unit("1"), unit("m"))
    end
end

"""
    pg_units(key::String)

Returns a PMDUnit representing the units of any attribute.
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

    if startswith(key, "cov_")
        subkeys = split(key[5:end], "__")
        unit0 = PARTICLEGROUP_UNITS[subkeys[1]]
        unit1 = PARTICLEGROUP_UNITS[subkeys[2]]
        return multiply_units(unit0, unit1)
    end

    # Fields
    if startswith(key, "electricField")
        return unit("V/m")
    end
    if startswith(key, "magneticField")
        return unit("T")
    end
    if startswith(key, "bunching_phase")
        return unit("rad")
    end
    if startswith(key, "bunching")
        return unit("1")
    end

    throw(ArgumentError("No known unit for: $key"))
end

"""
    parse_bunching_str(s::String)

Parse a string of the form to extract the wavelength:
    'bunching_1.23e-4'
    'bunching_1.23e-4_nm'
    'bunching_phase_1.23e-4'

Returns the wavelength in meters.
"""
function parse_bunching_str(s::String)
    @assert startswith(s, "bunching_")

    # Remove bunching and phase prefixes
    s = replace(s, "bunching_" => "")
    s = replace(s, "phase_" => "")

    x = split(s, "_")

    wavelength = parse(Float64, x[1])

    if length(x) == 1
        factor = 1.0
    elseif length(x) == 2
        unit = x[2]
        if unit == "m"
            factor = 1.0
        elseif unit == "mm"
            factor = 1e-3
        elseif unit in ["µm", "um"]
            factor = 1e-6
        elseif unit == "nm"
            factor = 1e-9
        else
            throw(ArgumentError("Unparsable unit: $unit"))
        end
    else
        throw(ArgumentError("Cannot parse $s"))
    end

    return wavelength * factor
end

# HDF5 tools
"""
    write_unit_h5(h5, u::PMDUnit)

Writes a PMDUnit to an h5 handle.
"""
function write_unit_h5(h5, u::PMDUnit)
    h5.attrs["unitSI"] = u.unitSI
    h5.attrs["unitDimension"] = u.unitDimension
    h5.attrs["unitSymbol"] = u.unitSymbol
end

"""
    read_unit_h5(h5)

Reads unit data from an h5 handle and returns a PMDUnit object.
"""
function read_unit_h5(h5)
    a = h5.attrs

    unitSI = a["unitSI"]
    unitDimension = Tuple(a["unitDimension"])
    unitSymbol = get(a, "unitSymbol", "unknown")

    return PMDUnit(unitSymbol=unitSymbol, unitSI=unitSI, unitDimension=unitDimension)
end

"""
    read_dataset_and_unit_h5(h5, expected_unit=nothing; convert=true)

Reads a dataset that has openPMD unit attributes.

expected_unit can be a PMDUnit object, or a known unit str. Examples: "kg", "J", "eV"

If expected_unit is given, will check that the units are compatible.

If convert, the data will be returned with the expected_units.

Returns a tuple: (array, PMDUnit)
"""
function read_dataset_and_unit_h5(h5, expected_unit=nothing; convert=true)
    # Read the unit that is there
    u = read_unit_h5(h5)

    # Simple case
    if isnothing(expected_unit)
        return Array(h5), u
    end

    if expected_unit isa String
        # Try to get unit
        expected_unit = unit(expected_unit)
    end

    # Check dimensions
    du = divide_units(u, expected_unit)

    @assert du.unitDimension == (0, 0, 0, 0, 0, 0, 0) "incompatible units"

    if convert
        fac = du.unitSI
        return fac * Array(h5), expected_unit
    else
        return Array(h5), u
    end
end

"""
    write_dataset_and_unit_h5(h5, name::String, data, unit=nothing)

Writes data and PMDUnit to h5[name].
"""
function write_dataset_and_unit_h5(h5, name::String, data, unit=nothing)
    h5[name] = data

    if !isnothing(unit)
        write_unit_h5(h5[name], unit)
    end
end
