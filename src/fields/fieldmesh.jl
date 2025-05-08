# Constants
const AXIS_LABELS_FROM_GEOMETRY = Dict(
    "rectangular" => ("x", "y", "z"),
    "cylindrical" => ("r", "theta", "z")
)

# Type definitions
mutable struct FieldMesh
    _data::Dict{String,Any}
    
    function FieldMesh(h5=nothing, data=nothing)
        if !isnothing(h5)
            # Allow filename
            if isa(h5, String)
                fname = expanduser(expandvars(h5))
                @assert isfile(fname) "File does not exist: $fname"
                
                h5file = h5open(fname, "r")
                fp = field_paths(h5file)
                @assert length(fp) == 1 "Number of field paths in $h5: $(length(fp))"
                data = load_field_data_h5(h5file[fp[1]])
                close(h5file)
            else
                data = load_field_data_h5(h5)
            end
        else
            data = load_field_data_dict(data)
        end
        
        # Internal data
        new(Dict("_data" => data))
    end
end

# Direct access to internal data
attrs(fm::FieldMesh) = fm._data["attrs"]
components(fm::FieldMesh) = fm._data["components"]
data(fm::FieldMesh) = fm._data

# Basic properties
shape(fm::FieldMesh) = Tuple(attrs(fm)["gridSize"])
geometry(fm::FieldMesh) = attrs(fm)["gridGeometry"]

# Scale property
scale(fm::FieldMesh) = attrs(fm)["fieldScale"]
function set_scale!(fm::FieldMesh, val)
    attrs(fm)["fieldScale"] = val
end

# Phase property
function phase(fm::FieldMesh)
    """
    Returns the complex argument `phi = -2*pi*RFphase`
    to multiply the oscillating field by.
    """
    return -attrs(fm)["RFphase"] * 2π
end

function set_phase!(fm::FieldMesh, val)
    """
    Complex argument in radians
    """
    attrs(fm)["RFphase"] = -val / (2π)
end

# Factor property
function factor(fm::FieldMesh)
    """
    factor to multiply fields by, possibly complex.
    `factor = scale * exp(i*phase)`
    """
    return real_if_close(scale(fm) * exp(im * phase(fm)))
end

# Axis labels
axis_labels(fm::FieldMesh) = AXIS_LABELS_FROM_GEOMETRY[geometry(fm)]

function axis_index(fm::FieldMesh, key)
    """
    Returns axis index for a named axis label key.
    """
    for (i, name) in enumerate(axis_labels(fm))
        if name == key
            return i
        end
    end
    error("Axis not found: $key")
end

# Coordinate vectors
function coord_vecs(fm::FieldMesh)
    """
    Uses gridSpacing, gridSize, and gridOriginOffset to return coordinate vectors.
    """
    mins = attrs(fm)["gridOriginOffset"]
    maxs = deltas(fm) .* (attrs(fm)["gridSize"] .- 1) .+ mins
    return [range(mins[i], maxs[i], length=attrs(fm)["gridSize"][i]) for i in 1:3]
end

function coord_vec(fm::FieldMesh, key)
    """
    Gets the coordinate vector from a named axis key.
    """
    i = axis_index(fm, key)
    mins = attrs(fm)["gridOriginOffset"]
    maxs = deltas(fm) .* (attrs(fm)["gridSize"] .- 1) .+ mins
    return range(mins[i], maxs[i], length=attrs(fm)["gridSize"][i])
end

# Meshgrid
function meshgrid(fm::FieldMesh)
    """
    Uses coordinate vectors to produce a standard meshgrid.
    """
    vecs = coord_vecs(fm)
    return collect.(Iterators.product(vecs...))
end

# Min/Max/Delta properties
mins(fm::FieldMesh) = attrs(fm)["gridOriginOffset"]
deltas(fm::FieldMesh) = attrs(fm)["gridSpacing"]
maxs(fm::FieldMesh) = deltas(fm) .* (attrs(fm)["gridSize"] .- 1) .+ mins(fm)

# Frequency property
function frequency(fm::FieldMesh)
    if is_static(fm)
        return 0
    else
        return attrs(fm)["harmonic"] * attrs(fm)["fundamentalFrequency"]
    end
end

# Logical properties
function is_pure_electric(fm::FieldMesh)
    """
    Returns true if there are no non-zero magnetic field components
    """
    klist = [key for key in keys(components(fm)) if !component_is_zero(fm, key)]
    return all(startswith(key, "electric") for key in klist)
end

function is_pure_magnetic(fm::FieldMesh)
    """
    Returns true if there are no non-zero electric field components
    """
    klist = [key for key in keys(components(fm)) if !component_is_zero(fm, key)]
    return all(startswith(key, "magnetic") for key in klist)
end

is_static(fm::FieldMesh) = attrs(fm)["harmonic"] == 0

function component_is_zero(fm::FieldMesh, key)
    """
    Returns true if all elements in a component are zero.
    """
    a = getindex(fm, key)
    return !any(a)
end

# Interpolation
function interpolator(fm::FieldMesh, key)
    """
    Returns an interpolator for a given field key.
    """
    field = getindex(fm, key)
    itp = interpolate(field, BSpline(Linear()))
    return extrapolate(itp, Flat())
end

function interpolate(fm::FieldMesh, key, points)
    """
    Interpolates the field data for the given key at specified points.
    """
    points = collect(points)
    
    # Convenience for a single point
    if ndims(points) == 1
        return interpolator(fm, key)([points])[1]
    end
    
    return interpolator(fm, key)(points)
end

# Helper functions
function real_if_close(x)
    """
    Returns real part if complex part is close to zero
    """
    if abs(imag(x)) < 1e-10
        return real(x)
    end
    return x
end

# Field component accessors
function scaled_component(fm::FieldMesh, key)
    """
    Returns a component scaled by the complex factor
    factor = scale*exp(i*phase)
    """
    if haskey(components(fm), key)
        dat = components(fm)[key]
    # Aliases
    elseif haskey(COMPONENT_FROM_ALIAS, key)
        comp = COMPONENT_FROM_ALIAS[key]
        if haskey(components(fm), comp)
            dat = components(fm)[comp]
        else
            # Component not present, make zeros
            return zeros(shape(fm))
        end
    else
        error("Component not available: $key")
    end

    # Multiply by scale factor
    factor_val = factor(fm)

    if factor_val != 1
        return factor_val * dat
    else
        return dat
    end
end

# Convenient properties
r(fm::FieldMesh) = coord_vec(fm, "r")
theta(fm::FieldMesh) = coord_vec(fm, "theta")
z(fm::FieldMesh) = coord_vec(fm, "z")

# Field components
function Bx(fm::FieldMesh)
    return scaled_component(fm, "Bx")
end

function By(fm::FieldMesh)
    return scaled_component(fm, "By")
end

function Bz(fm::FieldMesh)
    return scaled_component(fm, "Bz")
end

function Br(fm::FieldMesh)
    return scaled_component(fm, "Br")
end

function Btheta(fm::FieldMesh)
    return scaled_component(fm, "Btheta")
end

function Ex(fm::FieldMesh)
    return scaled_component(fm, "Ex")
end

function Ey(fm::FieldMesh)
    return scaled_component(fm, "Ey")
end

function Ez(fm::FieldMesh)
    return scaled_component(fm, "Ez")
end

function Er(fm::FieldMesh)
    return scaled_component(fm, "Er")
end

function Etheta(fm::FieldMesh)
    return scaled_component(fm, "Etheta")
end

function B(fm::FieldMesh)
    if geometry(fm) == "cylindrical"
        if is_static(fm)
            return hypot.(Br(fm), Bz(fm))
        else
            return abs.(Btheta(fm))
        end
    else
        error("Unknown geometry: $(geometry(fm))")
    end
end

function E(fm::FieldMesh)
    if geometry(fm) == "cylindrical"
        return hypot.(abs.(Er(fm)), abs.(Ez(fm)))
    else
        error("Unknown geometry: $(geometry(fm))")
    end
end

# Plotting
function plot(fm::FieldMesh; component=nothing, time=nothing, axes=nothing, cmap=nothing, return_figure=false, nice=true, kwargs...)
    if geometry(fm) == "cylindrical"
        return plot_fieldmesh_cylindrical_2d(
            fm,
            component=component,
            time=time,
            axes=axes,
            return_figure=return_figure,
            cmap=cmap,
            kwargs...
        )
    elseif geometry(fm) == "rectangular"
        return plot_fieldmesh_rectangular_2d(
            fm,
            component=component,
            time=time,
            axes=axes,
            return_figure=return_figure,
            nice=nice,
            cmap=cmap,
            kwargs...
        )
    else
        error("Geometry $(geometry(fm)) not implemented")
    end
end

function plot_onaxis(fm::FieldMesh, args...; kwargs...)
    if geometry(fm) == "cylindrical"
        return plot_fieldmesh_cylindrical_1d(fm, args...; kwargs...)
    elseif geometry(fm) == "rectangular"
        return plot_fieldmesh_rectangular_1d(fm, args...; kwargs...)
    else
        error("Unsupported geometry for plot_onaxis: $(geometry(fm))")
    end
end

# I/O functions
function write(fm::FieldMesh, h5, name=nothing)
    """
    Writes openPMD-beamphysics format to an open h5 handle, or new file if h5 is a str.
    """
    if isa(h5, String)
        fname = expanduser(expandvars(h5))
        h5file = h5open(fname, "w")
        pmd_field_init(h5file, externalFieldPath="/ExternalFieldPath/%T/")
        g = create_group(h5file, "/ExternalFieldPath/1/")
    else
        g = h5
    end

    write_pmd_field(g, data(fm), name=name)
end

# Utility functions
function get_operator(key)
    """
    Check if a key starts with `re_`, `im_`, `abs_`
    returns operator, newkey
    """
    # Check for simple operators
    if startswith(key, "re_")
        operator = real
        newkey = key[4:end]
    elseif startswith(key, "im_")
        operator = imag
        newkey = key[4:end]
    elseif startswith(key, "abs_")
        operator = abs
        newkey = key[5:end]
    else
        operator = nothing
        newkey = key
    end

    return operator, newkey
end

function load_field_data_h5(h5; verbose=true)
    """
    Loads field data from an HDF5 file.
    If `attrs['dataOrder'] == 'F'`, will transpose.
    If `attrs['harmonic'] == 0`, components will be cast to real by `real`
    """
    data = Dict("components" => Dict())

    # Load attributes
    attrs = decode_attrs(attributes(h5))
    data["attrs"] = attrs

    # Loop over records and components
    for (g, comps) in field_record_components
        if !haskey(h5, g)
            continue
        end

        # Extract axis labels for data transposing
        axis_labels = Tuple([decode_attr(a) for a in attributes(h5)["axisLabels"]])

        # Get the full openPMD unitDimension
        required_dim = expected_record_unit_dimension[g]

        # Filter out only the components that we have
        comps = [comp for comp in comps if haskey(h5[g], comp)]

        for comp in comps
            name = "$g/$comp"

            # Handle legacy format
            if is_legacy_fortran_data_ordering(attributes(h5[name]))
                if "r" in comps
                    axis_labels = ("z", "theta", "r")
                else
                    axis_labels = ("z", "y", "x")
                end
            end

            cdat = component_data(h5[name], axis_labels=axis_labels)

            # Check dimensions
            dim = attributes(h5[name])["unitDimension"]
            @assert all(dim .== required_dim) "$name with dimension $required_dim expected for $name, found: $dim"

            # Check shape
            s1 = Tuple(attrs["gridSize"])
            s2 = size(cdat)
            @assert s1 == s2 "Expected shape: $s1 != found shape: $s2"

            # Static fields should be real
            if attrs["harmonic"] == 0
                cdat = real(cdat)
            end

            # Finally set
            data["components"][name] = cdat
        end
    end

    return data
end

function load_field_data_dict(data_dict; verbose=true)
    """
    Similar to load_field_data_h5, but from a dict.
    This cannot do unit checking.
    """
    # The output dict
    data = Dict()

    # Load attributes
    attrs = decode_attrs(data_dict["attrs"])
    data["attrs"] = attrs

    # Go through components. Allow aliases
    comp = data["components"] = Dict()
    for (k, v) in data_dict["components"]
        if haskey(COMPONENT_ALIAS, k)
            comp[k] = v
        elseif haskey(COMPONENT_FROM_ALIAS, k)
            k = COMPONENT_FROM_ALIAS[k]
            @assert !haskey(data, k)
            comp[k] = v
        else
            error("Unallowed component: $k")
        end
    end

    return data
end

# Constructors
function from_ansys_ascii_3d(; efile=nothing, hfile=nothing, frequency=nothing)
    """
    Class method to return a FieldMesh from ANSYS ASCII files.
    """
    if isnothing(frequency)
        error("Please provide a frequency")
    end

    data = read_ansys_ascii_3d_fields(efile, hfile, frequency=frequency)
    return FieldMesh(data=data)
end

function from_cst_3d(field_file1, field_file2=nothing, frequency=0)
    if !isnothing(field_file2)
        # field_file1 -> efile, field_file2 -> hfile
        data = read_cst_ascii_3d_complex_fields(
            field_file1, field_file2, frequency=frequency, harmonic=1
        )
    else
        data = read_cst_ascii_3d_static_field(field_file1)
    end

    return FieldMesh(data=data)
end

function from_astra_3d(common_filename, frequency=0)
    """
    Class method to parse multiple 3D astra fieldmap files,
    based on the common filename.
    """
    data = read_astra_3d_fieldmaps(common_filename, frequency=frequency)
    return FieldMesh(data=data)
end

function from_impact_emfield_cartesian(filename, frequency=0, eleAnchorPt="beginning")
    """
    Class method to read an Impact-T style 1Tv3.T7 file corresponding to
    the `111: EMfldCart` element.
    """
    attrs, components = parse_impact_emfield_cartesian(filename)

    # These aren't in the file, they must be added
    attrs["fundamentalFrequency"] = frequency
    if frequency == 0
        attrs["harmonic"] = 0
    end

    attrs["eleAnchorPt"] = eleAnchorPt
    return FieldMesh(data=Dict("attrs" => attrs, "components" => components))
end

function from_superfish(filename; type=nothing, geometry="cylindrical")
    """
    Class method to parse a superfish T7 style file.
    """
    data = read_superfish_t7(filename, type=type, geometry=geometry)
    return FieldMesh(data=data)
end

function from_onaxis(; z=nothing, Bz=nothing, Ez=nothing, frequency=0, harmonic=nothing, eleAnchorPt="beginning")
    """
    Create a FieldMesh from on-axis field data.
    """
    # Get spacing
    nz = length(z)
    dz = diff(z)
    if !all(isapprox.(dz, dz[1]))
        error("Irregular spacing not implemented")
    end
    dz = dz[1]

    components = Dict()
    if !isnothing(Ez)
        Ez = squeeze(Ez)
        if ndims(Ez) != 1
            error("Ez ndim = $(ndims(Ez)) must be 1")
        end
        components["electricField/z"] = reshape(Ez, 1, 1, length(Ez))
    end

    if !isnothing(Bz)
        Bz = squeeze(Bz)
        if ndims(Bz) != 1
            error("Bz ndim = $(ndims(Bz)) must be 1")
        end
        components["magneticField/z"] = reshape(Bz, 1, 1, length(Bz))
    end

    if isnothing(Bz) && isnothing(Ez)
        error("Please enter Ez or Bz")
    end

    # Handle harmonic options
    if frequency == 0
        harmonic = 0
    elseif isnothing(harmonic)
        harmonic = 1
    end

    attrs = Dict(
        "eleAnchorPt" => eleAnchorPt,
        "gridGeometry" => "cylindrical",
        "axisLabels" => ["r", "theta", "z"],
        "gridLowerBound" => [0, 0, 0],
        "gridOriginOffset" => [0.0, 0.0, minimum(z)],
        "gridSpacing" => [0.0, 0.0, dz],
        "gridSize" => [1, 1, nz],
        "harmonic" => harmonic,
        "fundamentalFrequency" => frequency,
        "RFphase" => 0,
        "fieldScale" => 1.0
    )

    data = Dict("attrs" => attrs, "components" => components)
    return FieldMesh(data=data)
end

# Additional methods
function Base.getindex(fm::FieldMesh, key)
    """
    Returns component data from a key
    """
    if key in ["r", "theta", "z"]
        return coord_vec(fm, key)
    end

    # Raw components
    if haskey(components(fm), key)
        return components(fm)[key]
    end

    # Check for operators
    operator, key = get_operator(key)

    # Scaled components
    if key == "E"
        dat = E(fm)
    elseif key == "B"
        dat = B(fm)
    else
        dat = scaled_component(fm, key)
    end

    if !isnothing(operator)
        dat = operator(dat)
    end

    return dat
end

function Base.copy(fm::FieldMesh)
    """Returns a deep copy"""
    return deepcopy(fm)
end

function Base.show(io::IO, fm::FieldMesh)
    memloc = string(objectid(fm), base=16)
    print(io, "<FieldMesh with $(geometry(fm)) geometry and $(shape(fm)) shape at $memloc>")
end