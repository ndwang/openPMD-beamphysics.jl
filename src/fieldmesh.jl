"""
FieldMesh: openPMD-beamphysics External Field Mesh data structure.

Represents gridded electromagnetic field data in rectangular or cylindrical geometry,
following the openPMD-beamphysics standard.
"""

# ── Geometry helpers ──────────────────────────────────────────────────────

const AXIS_LABELS_FROM_GEOMETRY = Dict{String, NTuple{3, String}}(
    "rectangular" => ("x", "y", "z"),
    "cylindrical" => ("r", "theta", "z"),
)

# ── FieldMesh struct ──────────────────────────────────────────────────────

"""
    FieldMesh

Gridded electromagnetic field data in openPMD-beamphysics format.

# Fields
- `geometry`: `"rectangular"` or `"cylindrical"`
- `axis_labels`: e.g. `("r", "theta", "z")`
- `grid_origin`: physical coordinate of grid origin
- `grid_spacing`: spacing per axis
- `grid_size`: number of grid points per axis
- `ele_anchor_pt`: element anchor point (`"beginning"`, `"center"`, `"end"`)
- `harmonic`: 0 for static, >0 for RF
- `fundamental_frequency`: in Hz
- `rf_phase`: RF phase in units of 2π (openPMD `RFphase`)
- `field_scale`: overall field scale factor
- `components`: Dict mapping openPMD names (e.g. `"magneticField/z"`) to 3D arrays

# Constructors
- `FieldMesh(filename::String)` — from HDF5 file
- `FieldMesh(h5::HDF5.Group)` — from open HDF5 group
- Keyword constructor for programmatic creation

# Examples
```julia
fm = FieldMesh("solenoid.h5")
fm.geometry        # "cylindrical"
fm["Bz"]           # scaled Bz component
coord_vec(fm, "z") # z coordinate vector
```
"""
mutable struct FieldMesh
    geometry::String
    axis_labels::NTuple{3, String}
    grid_origin::NTuple{3, Float64}
    grid_spacing::NTuple{3, Float64}
    grid_size::NTuple{3, Int}
    ele_anchor_pt::String
    harmonic::Int
    fundamental_frequency::Float64
    rf_phase::Float64
    field_scale::Float64
    components::Dict{String, Array{Float64, 3}}
end

# ── Constructors ──────────────────────────────────────────────────────────

"""
    FieldMesh(filename::String)

Load a FieldMesh from an openPMD-beamphysics HDF5 file.
"""
function FieldMesh(filename::String)
    h5open(filename, "r") do h5
        fp = field_paths(h5)
        isempty(fp) && error("No field paths found in $filename")
        length(fp) > 1 && error("Multiple field paths in $filename: $fp. Use FieldMesh(h5group) to select one.")
        return FieldMesh(h5[fp[1]])
    end
end

"""
    FieldMesh(h5::HDF5.Group)

Load a FieldMesh from an open HDF5 group.
"""
function FieldMesh(h5::HDF5.Group)
    a, _ = load_field_attrs(attrs(h5))

    geometry = String(a["gridGeometry"])
    axis_labels = NTuple{3, String}(String.(a["axisLabels"]))
    grid_origin = NTuple{3, Float64}(a["gridOriginOffset"])
    grid_spacing = NTuple{3, Float64}(a["gridSpacing"])
    grid_size = NTuple{3, Int}(a["gridSize"])
    ele_anchor_pt = String(a["eleAnchorPt"])
    harmonic = _to_scalar(Int, a["harmonic"])
    fundamental_frequency = Float64(get(a, "fundamentalFrequency", 0))
    rf_phase = Float64(get(a, "RFphase", 0))
    field_scale = _to_scalar(Float64, a["fieldScale"])

    components = Dict{String, Array{Float64, 3}}()
    for (record, comp_names) in FIELD_RECORD_COMPONENTS
        haskey(h5, record) || continue
        for cname in comp_names
            haskey(h5[record], cname) || continue
            key = record * "/" * cname
            cdat = component_data(h5[key]; axis_labels=axis_labels)
            # Ensure 3D Float64 array
            cdat = Array{Float64, 3}(real(cdat))
            components[key] = cdat
        end
    end

    return FieldMesh(geometry, axis_labels, grid_origin, grid_spacing, grid_size,
                     ele_anchor_pt, harmonic, fundamental_frequency, rf_phase,
                     field_scale, components)
end

"""Extract a scalar from a value that might be a 1-element array."""
_to_scalar(T, x::AbstractArray) = T(x[1])
_to_scalar(T, x) = T(x)

# ── Grid accessors ────────────────────────────────────────────────────────

"""
    axis_index(fm::FieldMesh, key::String) -> Int

Return the axis index (1-based) for a named axis label.
"""
function axis_index(fm::FieldMesh, key::String)
    for i in 1:3
        fm.axis_labels[i] == key && return i
    end
    throw(ArgumentError("Axis not found: $key. Available: $(fm.axis_labels)"))
end

"""
    mins(fm::FieldMesh) -> NTuple{3, Float64}

Grid minimum coordinates.
"""
mins(fm::FieldMesh) = fm.grid_origin

"""
    deltas(fm::FieldMesh) -> NTuple{3, Float64}

Grid spacing.
"""
deltas(fm::FieldMesh) = fm.grid_spacing

"""
    maxs(fm::FieldMesh) -> NTuple{3, Float64}

Grid maximum coordinates.
"""
function maxs(fm::FieldMesh)
    return ntuple(i -> fm.grid_origin[i] + fm.grid_spacing[i] * (fm.grid_size[i] - 1), 3)
end

"""
    coord_vec(fm::FieldMesh, label::String) -> Vector{Float64}

Coordinate vector along the named axis.
"""
function coord_vec(fm::FieldMesh, label::String)
    i = axis_index(fm, label)
    return range(fm.grid_origin[i], step=fm.grid_spacing[i], length=fm.grid_size[i]) |> collect
end

"""
    coord_vecs(fm::FieldMesh) -> Tuple of Vector{Float64}

Coordinate vectors for all three axes.
"""
coord_vecs(fm::FieldMesh) = Tuple(coord_vec(fm, l) for l in fm.axis_labels)

# Per-axis delta accessors
for name in ("x", "y", "z", "r", "theta")
    fname = Symbol("d" * name)
    @eval begin
        """    $($fname)(fm::FieldMesh) — grid spacing in $($name)."""
        $fname(fm::FieldMesh) = fm.grid_spacing[axis_index(fm, $name)]
    end
end

# Per-axis min/max accessors
for name in ("x", "y", "z", "r", "theta")
    fname_min = Symbol(name * "min")
    fname_max = Symbol(name * "max")
    @eval begin
        """    $($fname_min)(fm::FieldMesh) — grid minimum in $($name)."""
        $fname_min(fm::FieldMesh) = mins(fm)[axis_index(fm, $name)]
        """    $($fname_max)(fm::FieldMesh) — grid maximum in $($name)."""
        $fname_max(fm::FieldMesh) = maxs(fm)[axis_index(fm, $name)]
    end
end

# ── Min/Max setters ───────────────────────────────────────────────────────

"""
    set_min!(fm::FieldMesh, label::String, val::Real)

Set the grid minimum along `label` by shifting the grid origin.
"""
function set_min!(fm::FieldMesh, label::String, val::Real)
    i = axis_index(fm, label)
    origin = collect(fm.grid_origin)
    origin[i] = Float64(val)
    fm.grid_origin = NTuple{3, Float64}(origin)
    return fm
end

"""
    set_max!(fm::FieldMesh, label::String, val::Real)

Set the grid maximum along `label` by shifting the grid origin.
"""
function set_max!(fm::FieldMesh, label::String, val::Real)
    i = axis_index(fm, label)
    current_max = maxs(fm)[i]
    origin = collect(fm.grid_origin)
    origin[i] += Float64(val) - current_max
    fm.grid_origin = NTuple{3, Float64}(origin)
    return fm
end

# ── Field scaling ─────────────────────────────────────────────────────────

"""
    scale(fm::FieldMesh) -> Float64

Field scale factor.
"""
scale(fm::FieldMesh) = fm.field_scale

"""
    phase(fm::FieldMesh) -> Float64

Complex phase argument `φ = -2π * RFphase` in radians.
"""
phase(fm::FieldMesh) = -fm.rf_phase * 2π

"""
    factor(fm::FieldMesh)

Complex factor to multiply fields by: `scale * exp(i * phase)`.
Returns real if phase is zero.
"""
function factor(fm::FieldMesh)
    p = phase(fm)
    f = fm.field_scale * exp(im * p)
    return p == 0 ? real(f) : f
end

"""
    frequency(fm::FieldMesh) -> Float64

Oscillation frequency in Hz. Returns 0 for static fields.
"""
frequency(fm::FieldMesh) = is_static(fm) ? 0.0 : fm.harmonic * fm.fundamental_frequency

# ── Boolean queries ───────────────────────────────────────────────────────

"""
    is_static(fm::FieldMesh) -> Bool

True if the field is static (harmonic == 0).
"""
is_static(fm::FieldMesh) = fm.harmonic == 0

"""
    is_pure_electric(fm::FieldMesh) -> Bool

True if there are no non-zero magnetic field components.
"""
function is_pure_electric(fm::FieldMesh)
    nonzero = [k for (k, v) in fm.components if any(!iszero, v)]
    return all(startswith(k, "electricField") for k in nonzero)
end

"""
    is_pure_magnetic(fm::FieldMesh) -> Bool

True if there are no non-zero electric field components.
"""
function is_pure_magnetic(fm::FieldMesh)
    nonzero = [k for (k, v) in fm.components if any(!iszero, v)]
    return all(startswith(k, "magneticField") for k in nonzero)
end

# ── Component access ──────────────────────────────────────────────────────

"""
    scaled_component(fm::FieldMesh, key::String) -> Array

Return a field component scaled by `factor(fm)`.
Accepts aliases like `"Bz"` or full names like `"magneticField/z"`.
"""
function scaled_component(fm::FieldMesh, key::String)
    comp_key = _resolve_component_key(fm, key)
    if comp_key === nothing
        return zeros(Float64, fm.grid_size)
    end
    dat = fm.components[comp_key]
    f = factor(fm)
    return f == 1 ? dat : f .* dat
end

"""Resolve an alias or full component key, returning nothing if not present."""
function _resolve_component_key(fm::FieldMesh, key::String)
    # Direct match
    haskey(fm.components, key) && return key
    # Alias lookup
    if haskey(COMPONENT_FROM_ALIAS, key)
        full = COMPONENT_FROM_ALIAS[key]
        haskey(fm.components, full) && return full
        return nothing
    end
    throw(ArgumentError("Unknown component: $key"))
end

"""
    _get_operator(key::String) -> (op, stripped_key)

Check if key starts with `re_`, `im_`, `abs_` and return the operator + stripped key.
"""
function _get_operator(key::String)
    startswith(key, "re_") && return (real, key[4:end])
    startswith(key, "im_") && return (imag, key[4:end])
    startswith(key, "abs_") && return (abs, key[5:end])
    return (nothing, key)
end

"""
    _field_magnitude(fm::FieldMesh, prefix::String) -> Array

Compute field magnitude for `"E"` or `"B"`.
"""
function _field_magnitude(fm::FieldMesh, prefix::String)
    if fm.geometry == "cylindrical"
        r_key = prefix == "E" ? "Er" : "Br"
        z_key = prefix == "E" ? "Ez" : "Bz"
        t_key = prefix == "E" ? "Etheta" : "Btheta"
        fr = scaled_component(fm, r_key)
        fz = scaled_component(fm, z_key)
        ft = scaled_component(fm, t_key)
        return @. sqrt(abs(fr)^2 + abs(fz)^2 + abs(ft)^2)
    elseif fm.geometry == "rectangular"
        x_key = prefix == "E" ? "Ex" : "Bx"
        y_key = prefix == "E" ? "Ey" : "By"
        z_key = prefix == "E" ? "Ez" : "Bz"
        fx = scaled_component(fm, x_key)
        fy = scaled_component(fm, y_key)
        fz = scaled_component(fm, z_key)
        return @. sqrt(abs(fx)^2 + abs(fy)^2 + abs(fz)^2)
    else
        throw(ArgumentError("Unknown geometry: $(fm.geometry)"))
    end
end

# ── getindex ──────────────────────────────────────────────────────────────

"""
    fm[key::String]

Access field data by key.

- `fm["Bz"]`, `fm["Ez"]` — scaled field components
- `fm["re_Bz"]`, `fm["im_Ez"]`, `fm["abs_Ez"]` — with operators
- `fm["r"]`, `fm["z"]` — coordinate vectors
- `fm["E"]`, `fm["B"]` — field magnitudes
"""
function Base.getindex(fm::FieldMesh, key::String)
    # Coordinate vectors
    if key in fm.axis_labels
        return coord_vec(fm, key)
    end

    # Check for operators
    op, stripped = _get_operator(key)

    # Field magnitudes
    dat = if stripped == "E"
        _field_magnitude(fm, "E")
    elseif stripped == "B"
        _field_magnitude(fm, "B")
    else
        scaled_component(fm, stripped)
    end

    return op === nothing ? dat : op.(dat)
end

# ── I/O ───────────────────────────────────────────────────────────────────

"""
    _to_data_dict(fm::FieldMesh) -> Dict

Convert FieldMesh to the data dict format expected by `write_pmd_field`.
"""
function _to_data_dict(fm::FieldMesh)
    attrs = Dict{String, Any}(
        "gridGeometry" => fm.geometry,
        "axisLabels" => collect(fm.axis_labels),
        "gridOriginOffset" => collect(fm.grid_origin),
        "gridSpacing" => collect(fm.grid_spacing),
        "gridSize" => collect(fm.grid_size),
        "gridLowerBound" => [0, 0, 0],
        "eleAnchorPt" => fm.ele_anchor_pt,
        "harmonic" => fm.harmonic,
        "fundamentalFrequency" => fm.fundamental_frequency,
        "RFphase" => fm.rf_phase,
        "fieldScale" => fm.field_scale,
    )
    components = Dict{String, Array}(k => v for (k, v) in fm.components)
    return Dict{String, Any}("attrs" => attrs, "components" => components)
end

"""
    write_fieldmesh(filename::String, fm::FieldMesh)

Write a FieldMesh to a new openPMD-beamphysics HDF5 file.
"""
function write_fieldmesh(filename::String, fm::FieldMesh)
    h5open(filename, "w") do h5
        pmd_field_init(h5)
        g = create_group(h5, "ExternalFieldPath/1")
        write_pmd_field(g, _to_data_dict(fm))
    end
end

"""
    write_fieldmesh(h5, fm::FieldMesh; name=nothing)

Write a FieldMesh to an open HDF5 group.
"""
function write_fieldmesh(h5, fm::FieldMesh; name=nothing)
    write_pmd_field(h5, _to_data_dict(fm); name=name)
end

# ── Base overloads ────────────────────────────────────────────────────────

function Base.show(io::IO, fm::FieldMesh)
    print(io, "FieldMesh($(fm.geometry), $(fm.grid_size))")
end

function Base.:(==)(a::FieldMesh, b::FieldMesh)
    a.geometry == b.geometry || return false
    a.axis_labels == b.axis_labels || return false
    a.grid_origin == b.grid_origin || return false
    a.grid_spacing == b.grid_spacing || return false
    a.grid_size == b.grid_size || return false
    a.ele_anchor_pt == b.ele_anchor_pt || return false
    a.harmonic == b.harmonic || return false
    a.fundamental_frequency == b.fundamental_frequency || return false
    a.rf_phase == b.rf_phase || return false
    a.field_scale == b.field_scale || return false
    keys(a.components) == keys(b.components) || return false
    return all(a.components[k] == b.components[k] for k in keys(a.components))
end

Base.copy(fm::FieldMesh) = FieldMesh(
    fm.geometry, fm.axis_labels, fm.grid_origin, fm.grid_spacing, fm.grid_size,
    fm.ele_anchor_pt, fm.harmonic, fm.fundamental_frequency, fm.rf_phase,
    fm.field_scale, Dict(k => copy(v) for (k, v) in fm.components)
)
