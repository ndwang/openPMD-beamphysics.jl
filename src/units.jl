"""
Self-contained units for the openPMD beamphysics records.

Replaces Unitful.jl with a lightweight PMDUnit struct that carries a symbol string,
an SI conversion factor, and the openPMD 7-tuple dimension vector (L,M,T,I,Θ,N,J).
"""

# ── PMDUnit struct ──────────────────────────────────────────────────────────

"""
    PMDUnit(symbol, scale, dims)

A physical unit described by its symbol string, SI scale factor, and openPMD
7-tuple dimension vector `(length, mass, time, current, temperature, amount, luminosity)`.

Dimensions are `Float64` to support fractional exponents (e.g. `m^(1/2)`).
"""
struct PMDUnit
    symbol::String
    scale::Float64
    dims::NTuple{7,Float64}

    function PMDUnit(symbol::String, scale::Float64, dims::NTuple{7,Float64})
        # Normalize -0.0 → +0.0 so that isequal / hash behave consistently
        new(symbol, scale, ntuple(i -> dims[i] + 0.0, 7))
    end
end

# Accept any Real for scale
PMDUnit(sym::String, scale::Real, dims::NTuple{7,Float64}) =
    PMDUnit(sym, Float64(scale), dims)

# Accept any iterable convertible to NTuple{7,Float64}
PMDUnit(sym::String, scale::Real, dims) =
    PMDUnit(sym, Float64(scale), NTuple{7,Float64}(dims))

# ── Dimension dictionaries ──────────────────────────────────────────────────

const DIMENSION = Dict{String, NTuple{7,Float64}}(
    "1"                  => (0., 0., 0., 0., 0., 0., 0.),
    "length"             => (1., 0., 0., 0., 0., 0., 0.),
    "mass"               => (0., 1., 0., 0., 0., 0., 0.),
    "time"               => (0., 0., 1., 0., 0., 0., 0.),
    "current"            => (0., 0., 0., 1., 0., 0., 0.),
    "temperature"        => (0., 0., 0., 0., 1., 0., 0.),
    "mol"                => (0., 0., 0., 0., 0., 1., 0.),
    "luminous"           => (0., 0., 0., 0., 0., 0., 1.),
    "charge"             => (0., 0., 1., 1., 0., 0., 0.),
    "electric_field"     => (1., 1.,-3.,-1., 0., 0., 0.),
    "electric_potential" => (2., 1.,-3.,-1., 0., 0., 0.),
    "magnetic_field"     => (0., 1.,-2.,-1., 0., 0., 0.),
    "velocity"           => (1., 0.,-1., 0., 0., 0., 0.),
    "energy"             => (2., 1.,-2., 0., 0., 0., 0.),
    "momentum"           => (1., 1.,-1., 0., 0., 0., 0.),
)

const DIMENSION_NAME = Dict{NTuple{7,Float64}, String}(v => k for (k, v) in DIMENSION)

const SI_SYMBOL = Dict{String, String}(
    "1"                  => "1",
    "length"             => "m",
    "mass"               => "kg",
    "time"               => "s",
    "current"            => "A",
    "temperature"        => "K",
    "mol"                => "mol",
    "luminous"           => "cd",
    "charge"             => "C",
    "electric_field"     => "V/m",
    "electric_potential"  => "V",
    "velocity"           => "m/s",
    "energy"             => "J",
    "momentum"           => "kg⋅m/s",
    "magnetic_field"     => "T",
)

const SI_NAME = Dict{String, String}(v => k for (k, v) in SI_SYMBOL)

# ── 3-arg convenience constructor (with dimension name) ─────────────────────

"""
    PMDUnit(symbol, scale, dim_name::String)

Construct a `PMDUnit` by looking up `dim_name` in the `DIMENSION` table.
"""
PMDUnit(sym::String, scale::Real, dim_name::String) =
    PMDUnit(sym, Float64(scale), DIMENSION[dim_name])

# ── dimension / dimension_name ──────────────────────────────────────────────

"""
    dimension(name::String) -> NTuple{7,Float64}

Get the dimension tuple for a named quantity.
"""
function dimension(name::String)
    haskey(DIMENSION, name) || error(
        "Invalid unit dimension string: $name. " *
        "Valid options are: $(join(keys(DIMENSION), ", "))")
    return DIMENSION[name]
end

"""
    dimension_name(dims::NTuple{7,Float64}) -> String

Get the human-readable name for a dimension tuple.
"""
function dimension_name(dims::NTuple{7,Float64})
    haskey(DIMENSION_NAME, dims) || error("Unknown dimension tuple: $dims")
    return DIMENSION_NAME[dims]
end

"""
    dimension_name(u::PMDUnit) -> String
"""
dimension_name(u::PMDUnit) = dimension_name(u.dims)

# ── Predicates ──────────────────────────────────────────────────────────────

is_dimensionless(u::PMDUnit) = all(d == 0.0 for d in u.dims)
is_identity(u::PMDUnit)      = is_dimensionless(u) && u.scale == 1.0

# ── Base overloads ──────────────────────────────────────────────────────────

Base.show(io::IO, u::PMDUnit) = print(io, u.symbol)

Base.:(==)(a::PMDUnit, b::PMDUnit) =
    a.symbol == b.symbol && a.scale == b.scale && a.dims == b.dims

Base.hash(u::PMDUnit, h::UInt) =
    hash(u.symbol, hash(u.scale, hash(u.dims, h)))

# ── Arithmetic ──────────────────────────────────────────────────────────────

function _format_power(p::Real)
    if isinteger(p) && abs(p) < 100
        return string(Int(p))
    elseif p == 0.5
        return "(1/2)"
    elseif p == -0.5
        return "(-1/2)"
    elseif p == 1.5
        return "(3/2)"
    elseif p == -1.5
        return "(-3/2)"
    else
        return string(p)
    end
end

"""
    power_unit(u::PMDUnit, p::Real) -> PMDUnit

Raise unit to power `p`.
"""
function power_unit(u::PMDUnit, p::Real)
    p == 0 && return PMDUnit("", 1.0, ntuple(_ -> 0.0, 7))
    p == 1 && return u

    if isempty(u.symbol)
        return PMDUnit("", u.scale^p, NTuple{7,Float64}(p .* u.dims))
    end

    ps = _format_power(p)
    sym = occursin(r"[*/]", u.symbol) ? "($(u.symbol))^$ps" : "$(u.symbol)^$ps"
    return PMDUnit(sym, u.scale^p, NTuple{7,Float64}(p .* u.dims))
end

"""
    sqrt_unit(u::PMDUnit) -> PMDUnit

Take the square root of a unit (power 1/2).
"""
sqrt_unit(u::PMDUnit) = power_unit(u, 0.5)

function Base.:*(a::PMDUnit, b::PMDUnit)
    is_identity(a) && return b
    is_identity(b) && return a
    sym = a.symbol == b.symbol ? "($(a.symbol))^2" : "$(a.symbol)*$(b.symbol)"
    return PMDUnit(sym, a.scale * b.scale, NTuple{7,Float64}(a.dims .+ b.dims))
end

function Base.:/(a::PMDUnit, b::PMDUnit)
    if a.symbol == b.symbol && a.scale == b.scale
        return PMDUnit("", 1.0, ntuple(_ -> 0.0, 7))
    end
    is_identity(b) && return a
    sym = "$(a.symbol)/$(b.symbol)"
    return PMDUnit(sym, a.scale / b.scale, NTuple{7,Float64}(a.dims .- b.dims))
end

# ── Named-unit registry ────────────────────────────────────────────────────

const NAMED_UNITS = [
    PMDUnit("",       1.0,                "1"),
    PMDUnit("rad",    1.0,                "1"),
    PMDUnit("degree", π/180,              "1"),
    PMDUnit("m",      1.0,                "length"),
    PMDUnit("kg",     1.0,                "mass"),
    PMDUnit("g",      1e-3,               "mass"),
    PMDUnit("s",      1.0,                "time"),
    PMDUnit("A",      1.0,                "current"),
    PMDUnit("K",      1.0,                "temperature"),
    PMDUnit("mol",    1.0,                "mol"),
    PMDUnit("cd",     1.0,                "luminous"),
    PMDUnit("C",      1.0,                "charge"),
    PMDUnit("V/m",    1.0,                "electric_field"),
    PMDUnit("V",      1.0,                "electric_potential"),
    PMDUnit("c",      C_LIGHT,            "velocity"),
    PMDUnit("m/s",    1.0,                "velocity"),
    PMDUnit("eV",     E_CHARGE,           "energy"),
    PMDUnit("J",      1.0,                "energy"),
    PMDUnit("eV/c",   E_CHARGE/C_LIGHT,   "momentum"),
    PMDUnit("W",      1.0, (2., 1., -3., 0., 0., 0., 0.)),          # power
    PMDUnit("T",      1.0,                "magnetic_field"),
]

const KNOWN_UNIT = Dict{String,PMDUnit}(u.symbol => u for u in NAMED_UNITS)
KNOWN_UNIT["1"]       = KNOWN_UNIT[""]      # alias: "1" → dimensionless
KNOWN_UNIT["c_light"] = KNOWN_UNIT["c"]     # legacy alias

# ── Symbol parser ───────────────────────────────────────────────────────────

"""
    PMDUnit(symbol::String)

Construct a `PMDUnit` by parsing a unit symbol string such as `"eV/c"`,
`"m^(1/2)"`, or `"kg*m/s"`.
"""
function PMDUnit(symbol::String)
    haskey(KNOWN_UNIT, symbol) && return KNOWN_UNIT[symbol]
    return _parse_unit_symbol(symbol)
end

function _parse_unit_symbol(s::String)
    # Handle sqrt(…) wrapper
    m = match(r"^sqrt\((.+)\)$", s)
    if m !== nothing
        inner = PMDUnit(String(m[1]))
        r = sqrt_unit(inner)
        return PMDUnit(s, r.scale, r.dims)
    end

    scale = 1.0
    dims  = ntuple(_ -> 0.0, 7)
    sign  = 1          # +1 = multiply, -1 = divide
    token_start  = 1
    paren_depth  = 0

    for i in eachindex(s)
        c = s[i]
        if c == '('
            paren_depth += 1
        elseif c == ')'
            paren_depth -= 1
        elseif paren_depth == 0 && (c == '*' || c == '/')
            tok = strip(s[token_start:prevind(s, i)])
            if !isempty(tok)
                u = _parse_single_unit(String(tok))
                if sign == 1
                    scale *= u.scale
                    dims   = dims .+ u.dims
                else
                    scale /= u.scale
                    dims   = dims .- u.dims
                end
            end
            sign = c == '*' ? 1 : -1
            token_start = nextind(s, i)
        end
    end

    # Last token
    tok = strip(s[token_start:end])
    if !isempty(tok)
        u = _parse_single_unit(String(tok))
        if sign == 1
            scale *= u.scale
            dims   = dims .+ u.dims
        else
            scale /= u.scale
            dims   = dims .- u.dims
        end
    end

    return PMDUnit(s, scale, NTuple{7,Float64}(dims))
end

function _parse_single_unit(token::String)
    # Handle power notation: base^exp
    m = match(r"^(.+)\^(.+)$", token)
    if m !== nothing
        base  = String(m[1])
        power = _parse_power(String(m[2]))
        return power_unit(_lookup_unit(base), power)
    end
    return _lookup_unit(token)
end

function _parse_power(s::String)
    # (n/d) fractional form, e.g. "(1/2)", "(-3/2)"
    m = match(r"^\((-?\d+)/(\d+)\)$", s)
    m !== nothing && return parse(Float64, m[1]) / parse(Float64, m[2])
    return parse(Float64, s)
end

function _lookup_unit(name::String)
    haskey(KNOWN_UNIT, name) || error("Unknown unit: \"$name\"")
    return KNOWN_UNIT[name]
end

# ── SI prefix tables ───────────────────────────────────────────────────────

const PREFIX_FACTOR = Dict{String, Float64}(
    "yocto-" => 1e-24, "zepto-" => 1e-21, "atto-"  => 1e-18,
    "femto-" => 1e-15, "pico-"  => 1e-12, "nano-"  => 1e-9,
    "micro-" => 1e-6,  "milli-" => 1e-3,  "centi-" => 1e-2,
    "deci-"  => 1e-1,  "deca-"  => 1e1,   "hecto-" => 1e2,
    "kilo-"  => 1e3,   "mega-"  => 1e6,   "giga-"  => 1e9,
    "tera-"  => 1e12,  "peta-"  => 1e15,  "exa-"   => 1e18,
    "zetta-" => 1e21,  "yotta-" => 1e24,
)

const PREFIX = Dict{Float64, String}(v => k for (k, v) in PREFIX_FACTOR)

const SHORT_PREFIX_FACTOR = Dict{String, Float64}(
    "y"  => 1e-24, "z"  => 1e-21, "a"  => 1e-18,
    "f"  => 1e-15, "p"  => 1e-12, "n"  => 1e-9,
    "µ"  => 1e-6,  "m"  => 1e-3,  "c"  => 1e-2,
    "d"  => 1e-1,  ""   => 1,     "da" => 1e1,
    "h"  => 1e2,   "k"  => 1e3,   "M"  => 1e6,
    "G"  => 1e9,   "T"  => 1e12,  "P"  => 1e15,
    "E"  => 1e18,  "Z"  => 1e21,  "Y"  => 1e24,
)

const SHORT_PREFIX = Dict{Float64, String}(v => k for (k, v) in SHORT_PREFIX_FACTOR)

# ── nice_scale_prefix / nice_array / limits / plottable_array ───────────────

"""
    nice_scale_prefix(scale::Real) -> Tuple{Float64, String}

Return a nice scaling factor and the corresponding SI prefix string.
"""
function nice_scale_prefix(scale::Real)
    if scale == 0
        return 1.0, ""
    end

    p10 = log10(abs(scale))

    if p10 < -24
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

Scale an array and return `(scaled, factor, prefix)`.
"""
function nice_array(a::AbstractArray)
    if length(a) == 1
        x = a[1]
    else
        x = maximum((maximum(a) - minimum(a), abs(mean(a))))
    end
    fac, prefix = nice_scale_prefix(x)
    return a ./ fac, fac, prefix
end

"""
    nice_array(a) -> Tuple

Scalar version.
"""
function nice_array(a)
    fac, prefix = nice_scale_prefix(a)
    return a / fac, fac, prefix
end

"""
    limits(x::AbstractArray) -> Tuple
"""
limits(x::AbstractArray) = (minimum(x), maximum(x))

"""
    plottable_array(x; lim=nothing) -> Tuple

Scale an array for plotting, returning `(scaled, factor, prefix, xmin, xmax)`.
"""
function plottable_array(x::AbstractArray; lim::Union{Nothing, Tuple{Real, Real}}=nothing)
    xmin, xmax = lim === nothing ? limits(x) : lim
    _, factor, p1 = nice_array([xmin, xmax])
    return x ./ factor, factor, p1, xmin, xmax
end

# ── HDF5 I/O ───────────────────────────────────────────────────────────────

"""
    write_unit_h5(h5, u::PMDUnit)

Write a PMDUnit to an HDF5 handle as openPMD unit attributes.
"""
function write_unit_h5(h5, u::PMDUnit)
    attrs(h5)["unitSI"]        = u.scale
    attrs(h5)["unitDimension"] = collect(u.dims)
    attrs(h5)["unitSymbol"]    = u.symbol
end

"""
    read_unit_h5(h5) -> PMDUnit

Read unit metadata from an HDF5 handle and return a `PMDUnit`.
"""
function read_unit_h5(h5)
    sym = String(attrs(h5)["unitSymbol"])
    haskey(KNOWN_UNIT, sym) && return KNOWN_UNIT[sym]
    # Fallback: reconstruct from raw attributes
    return PMDUnit(sym, Float64(attrs(h5)["unitSI"]),
                   NTuple{7,Float64}(attrs(h5)["unitDimension"]))
end

"""
    read_dataset_and_unit_h5(h5)

Read a dataset that has openPMD unit attributes.
"""
function read_dataset_and_unit_h5(h5)
    u = read_unit_h5(h5)
    return Array(h5), u
end

"""
    read_dataset_and_unit_h5(h5, expected_unit::String; convert=true)

Read a dataset and optionally convert to `expected_unit`.
"""
function read_dataset_and_unit_h5(h5, expected_unit::String; convert::Bool=true)
    u  = read_unit_h5(h5)
    eu = PMDUnit(expected_unit)
    if u.dims != eu.dims
        error("Unit dimensions do not match: $(u.symbol) vs $(eu.symbol)")
    end
    if convert
        fac = u.scale / eu.scale
        return fac * Array(h5), eu
    else
        return Array(h5), u
    end
end

"""
    read_dataset_and_unit_h5(h5, expected_unit::PMDUnit; convert=true)

Read a dataset and optionally convert to `expected_unit`.
"""
function read_dataset_and_unit_h5(h5, expected_unit::PMDUnit; convert::Bool=true)
    u = read_unit_h5(h5)
    if u.dims != expected_unit.dims
        error("Unit dimensions do not match: $(u.symbol) vs $(expected_unit.symbol)")
    end
    if convert
        fac = u.scale / expected_unit.scale
        return fac * Array(h5), expected_unit
    else
        return Array(h5), u
    end
end

"""
    write_dataset_and_unit_h5(h5, name, data; unit=nothing)

Write data and unit to an HDF5 dataset.
"""
function write_dataset_and_unit_h5(h5, name::String, data; unit=nothing)
    h5[name] = data
    if unit !== nothing
        write_unit_h5(h5[name], unit)
    end
end

# ── ParticleGroup units ─────────────────────────────────────────────────────

const _DL = KNOWN_UNIT[""]       # dimensionless
const _M  = KNOWN_UNIT["m"]
const _S  = KNOWN_UNIT["s"]
const _EV = KNOWN_UNIT["eV"]
const _PC = KNOWN_UNIT["eV/c"]
const _RAD = KNOWN_UNIT["rad"]
const _COUL = KNOWN_UNIT["C"]
const _AMP  = KNOWN_UNIT["A"]
const _VM  = KNOWN_UNIT["V/m"]
const _T   = KNOWN_UNIT["T"]

const PARTICLEGROUP_UNITS = Dict{String, PMDUnit}(
    # Basic quantities
    "n_particle" => _DL,
    "status"     => _DL,
    "id"         => _DL,
    "n_alive"    => _DL,
    "n_dead"     => _DL,

    # Time
    "t"   => _S,
    "z/c" => _S,

    # Energy and mass
    "energy"                    => _EV,
    "kinetic_energy"            => _EV,
    "mass"                      => _EV,
    "higher_order_energy_spread"=> _EV,
    "higher_order_energy"       => _EV,

    # Momentum
    "px"     => _PC,
    "py"     => _PC,
    "pz"     => _PC,
    "p"      => _PC,
    "pr"     => _PC,
    "ptheta" => _PC,

    # Position
    "x"  => _M,
    "y"  => _M,
    "z"  => _M,
    "r"  => _M,
    "Jx" => _M,
    "Jy" => _M,

    # Dimensionless
    "beta"    => _DL,
    "beta_x"  => _DL,
    "beta_y"  => _DL,
    "beta_z"  => _DL,
    "gamma"   => _DL,
    "bunching" => _DL,

    # Angles
    "theta"          => _RAD,
    "bunching_phase" => _RAD,

    # Charge
    "charge"          => _COUL,
    "species_charge"  => _COUL,
    "weight"          => _COUL,
    "average_current" => _AMP,

    # Emittance
    "norm_emit_x"  => _M,
    "norm_emit_y"  => _M,
    "norm_emit_4d" => power_unit(_M, 2),

    # Angular momentum
    "Lz" => _M * _PC,

    # Angles
    "xp" => _RAD,
    "yp" => _RAD,

    # Normalized coordinates
    "x_bar"  => sqrt_unit(_M),
    "px_bar" => sqrt_unit(_M),
    "y_bar"  => sqrt_unit(_M),
    "py_bar" => sqrt_unit(_M),
)

# Field components
for component in ["", "x", "y", "z", "theta", "r"]
    PARTICLEGROUP_UNITS["E$(component)"] = _VM
    PARTICLEGROUP_UNITS["B$(component)"] = _T
end

# Twiss parameters
for plane in ("x", "y")
    PARTICLEGROUP_UNITS["twiss_alpha_$(plane)"]     = _DL
    PARTICLEGROUP_UNITS["twiss_etap_$(plane)"]      = _DL
    PARTICLEGROUP_UNITS["twiss_beta_$(plane)"]      = _M
    PARTICLEGROUP_UNITS["twiss_eta_$(plane)"]       = _M
    PARTICLEGROUP_UNITS["twiss_emit_$(plane)"]      = _M
    PARTICLEGROUP_UNITS["twiss_norm_emit_$(plane)"] = _M
    PARTICLEGROUP_UNITS["twiss_gamma_$(plane)"]     = power_unit(_M, -1)
end

"""
    pg_units(key::String) -> PMDUnit

Return the unit for any attribute in a ParticleGroup.
"""
function pg_units(key::String)
    haskey(PARTICLEGROUP_UNITS, key) && return PARTICLEGROUP_UNITS[key]

    # Operators: sigma_, mean_, min_, max_, ptp_, delta_
    for prefix in ["sigma_", "mean_", "min_", "max_", "ptp_", "delta_"]
        if startswith(key, prefix)
            return pg_units(key[length(prefix)+1:end])
        end
    end

    # Covariance
    if startswith(key, "cov_")
        subkeys = split(key[5:end], "__")
        length(subkeys) == 2 || error("Invalid covariance key format: $key")
        return PARTICLEGROUP_UNITS[subkeys[1]] * PARTICLEGROUP_UNITS[subkeys[2]]
    end

    # Fields
    startswith(key, "electricField") && return _VM
    startswith(key, "magneticField") && return _T
    startswith(key, "bunching_phase") && return _RAD
    startswith(key, "bunching") && return _DL

    error("No known unit for: $key")
end
