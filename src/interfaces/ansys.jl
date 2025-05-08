"""
    get_vec(x)

Helper function to get vector properties from a sorted array.
Returns min, max, spacing, and number of points.
"""
function get_vec(x)
    sx = unique(x)
    nx = length(sx)
    xlist = sort(sx)
    dx = diff(xlist)
    all(dx .≈ dx[1]) || error("Non-uniform spacing")
    dx = dx[1]
    return minimum(x), maximum(x), dx, nx
end

"""
    parse_ansys_ascii_3d(filePath, n_header=2)

Parses a single 3d field file.

The format is:
header1
header2
x y z re_fx im_fx re_fy im_fy re_fz im_fz
...

in C order

Ansys fields oscillate as:
    exp(i*omega*t)
Which is the opposite of openPMD-beamphysics's convention:
    exp(-i*omega*t)
"""
function parse_ansys_ascii_3d(filePath::String, n_header::Int=2)
    # Read data, skipping header
    dat = readdlm(filePath, skipstart=n_header)
    X = dat[:, 1]
    Y = dat[:, 2]
    Z = dat[:, 3]

    xmin, xmax, dx, nx = get_vec(X)
    ymin, ymax, dy, ny = get_vec(Y)
    zmin, zmax, dz, nz = get_vec(Z)

    shape = (nx, ny, nz)

    # - sign to convert to exp(-i omega t)
    Fx = reshape(complex.(dat[:, 4], -dat[:, 5]), shape)
    Fy = reshape(complex.(dat[:, 6], -dat[:, 7]), shape)
    Fz = reshape(complex.(dat[:, 8], -dat[:, 9]), shape)

    attrs = Dict{String,Any}(
        "gridOriginOffset" => (xmin, ymin, zmin),
        "gridSpacing" => (dx, dy, dz),
        "gridSize" => (nx, ny, nz)
    )

    components = Dict{String,Any}(
        "x" => Fx,
        "y" => Fy,
        "z" => Fz
    )

    return attrs, components
end

"""
    read_ansys_ascii_3d_fields(efile, hfile, frequency=0)

Reads electric and magnetic field data from ANSYS ASCII files.
"""
function read_ansys_ascii_3d_fields(efile::String, hfile::String, frequency::Float64=0.0)
    attrs1, components1 = parse_ansys_ascii_3d(efile)
    attrs2, components2 = parse_ansys_ascii_3d(hfile)

    # Check consistency
    for k in keys(attrs1)
        v1 = attrs1[k]
        v2 = attrs2[k]
        v1 == v2 || error("Inconsistent values for $k: $v1, $v2")
    end

    components = Dict{String,Any}()
    for k in keys(components1)
        components["electricField/$k"] = components1[k]
    end
    for k in keys(components2)
        components["magneticField/$k"] = components2[k] .* mu_0
    end

    attrs = attrs1
    attrs["eleAnchorPt"] = "beginning"
    attrs["gridGeometry"] = "rectangular"
    attrs["axisLabels"] = ("x", "y", "z")
    attrs["gridLowerBound"] = (0, 0, 0)
    attrs["harmonic"] = 1
    attrs["fundamentalFrequency"] = frequency

    return Dict{String,Any}("attrs" => attrs, "components" => components)
end 