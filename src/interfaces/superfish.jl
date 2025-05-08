"""
    write_fish_t7(fm, filePath, fmt="%10.8e", verbose=false)

Writes a T7 file from FISH t7data dict.

Input:
    fm: FieldMesh object
    filePath: requested filePath to write

See:
    superfish.parsers.parse_fish_t7
"""
function write_fish_t7(fm, filePath::String, fmt::String="%10.8e", verbose::Bool=false)
    fm.geometry == "cylindrical" || error("TODO: cartesian.")
    fm.frequency != 0 || error("Frequency must be non-zero.")

    rmin, _, zmin = fm.mins
    rmax, _, zmax = fm.maxs
    nr, _, nz = fm.shape

    # Collect these. Units are cm, MHz
    xmin = zmin * 100
    xmax = zmax * 100
    nx = nz
    ymin = rmin * 100
    ymax = rmax * 100
    ny = nr
    freq = fm.frequency * 1e-6

    # Get complex fields (helper function)
    Er, Ez, Btheta, _ = fish_complex_to_real_fields(fm, verbose=verbose)

    # Scale to Superfish units
    Er .*= 1e-6  # V/m -> MV/m
    Ez .*= 1e-6
    Hphi = Btheta ./ mu_0  # convert to H field, and in Superfish phase convention.
    E = hypot.(Er, Ez)

    # Write T7 ASCII
    header = """$xmin $xmax $(nx-1)
$freq
$ymin $ymax $(ny-1)"""

    # Unroll the arrays
    dat = hcat(reshape(Ez, nx * ny), reshape(Er, nx * ny), reshape(E, nx * ny), reshape(Hphi, nx * ny))

    open(filePath, "w") do io
        println(io, header)
        writedlm(io, dat, " ", fmt=fmt)
    end

    verbose && println("Superfish T7 file '$filePath' written for Fish problem.")

    return filePath
end

"""
    fish_complex_to_real_fields(fm, verbose=false)

Internal routine for .interfaces.superfish.write_fish_t7

Takes complex Ez, Er, Btheta and determines the complex angle to rotate Ez onto the real axis.

Returns rotated real fields:
    real(Er), real(Ez), -imag(Btheta), rot_angle
"""
function fish_complex_to_real_fields(fm, verbose::Bool=false)
    # Get complex fields
    Ez = fm.Ez[:, 1, :]
    Er = fm.Er[:, 1, :]
    Btheta = fm.Btheta[:, 1, :]

    # Get imaginary argument (angle)
    iangle = unique(mod.(angle.(Ez), π))

    # Make sure there is only one
    length(iangle) == 1 || error("Multiple angles found: $iangle")
    iangle = iangle[1]

    if iangle != 0
        rot = exp(-im * iangle)
        verbose && println("Rotating complex field by $(-iangle)")
    else
        rot = 1
    end

    Er = real.(Er .* rot)
    Ez = real.(Ez .* rot)
    Btheta = -imag.(Btheta .* rot)

    return Er, Ez, Btheta, -iangle
end

"""
    write_poisson_t7(fm, filePath, fmt="%10.8e", verbose=false)

Writes a T7 file from POISSON t7data dict.

See:
    superfish.parsers.parse_poisson_t7
"""
function write_poisson_t7(fm, filePath::String, fmt::String="%10.8e", verbose::Bool=false)
    fm.geometry == "cylindrical" || error("TODO: cartesian.")
    fm.is_static || error("Static fields are required for Poisson T7")

    rmin, _, zmin = fm.mins
    rmax, _, zmax = fm.maxs
    nr, _, nz = fm.shape

    # Collect these. Units are cm
    # Note: different from FISH!
    ymin = zmin * 100
    ymax = zmax * 100
    ny = nz
    xmin = rmin * 100
    xmax = rmax * 100
    nx = nr

    # Write T7 ASCII
    header = """$xmin $xmax $(nx-1)
$ymin $ymax $(ny-1)"""

    if fm.is_pure_electric
        kr = "Er"
        kz = "Ez"
        ftype = "electric"
        factor = 1
    else
        kr = "Br"
        kz = "Bz"
        ftype = "magnetic"
        factor = 1e4  # T->G
    end
    fr = real.(fm[kr][:, 1, :]) .* factor
    fz = real.(fm[kz][:, 1, :]) .* factor

    # Unroll the arrays
    dat = hcat(reshape(fr, nx * ny), reshape(fz, nx * ny))

    open(filePath, "w") do io
        println(io, header)
        writedlm(io, dat, " ", fmt=fmt)
    end

    verbose && println("Superfish T7 file '$filePath' written for $ftype Poisson problem.")

    return filePath
end

"""
    write_superfish_t7(fm, filePath, fmt="%10.8e", verbose=false)

Writes a Superfish T7 file. This is a simple wrapper for:
    write_fish_t7
    write_poisson_t7

If .is_static, a Poisson file is written. Otherwise a Fish file is written.

Parameters
----------
filePath: String
    File to write to

fmt: String, default = "%10.8e"
    Format to write numbers

Returns
-------
filePath: String
    File written (same as input)
"""
function write_superfish_t7(fm, filePath::String, fmt::String="%10.8e", verbose::Bool=false)
    if fm.is_static
        return write_poisson_t7(fm, filePath, fmt=fmt, verbose=verbose)
    else
        return write_fish_t7(fm, filePath, fmt=fmt, verbose=verbose)
    end
end

"""
    read_superfish_t7(filename, type=nothing, geometry="cylindrical")

Parses a T7 file written by Posson/Superfish.

Fish or Poisson T7 are automatically detected according to the second line.

For Poisson problems, the type must be specified.

Superfish fields oscillate as:
    Er, Ez ~ cos(wt)
    Hphi   ~ -sin(wt)

For complex fields oscillating as e^-iwt

    Re(Ex*e^-iwt)   ~ cos
    Re(-iB*e^-iwt) ~ -sin
and therefore B = -i * mu_0 * H_phi is the complex magnetic field in Tesla

Parameters:
----------
filename: String
    T7 filename to read
type: String, optional
    For Poisson files, required to be 'electric' or 'magnetic'.
    Not used for Fish files
geometry: String, optional
    field geometry, currently required to be the default: 'cylindrical'

Returns:
-------
fieldmesh_data: Dict
    Dictionary containing:
        attrs: Dict
        components: Dict

A FieldMesh object is instantiated from this as:
    FieldMesh(data=fieldmesh_data)
"""
function read_superfish_t7(filename::String, type=nothing, geometry::String="cylindrical")
    # ASCII parsing

    # Read header and data.
    # zmin(cm), zmax(cm), nx-1
    # freq(MHz)
    # ymin(cm), ymax(cm), ny-1
    open(filename, "r") do io
        line1 = split(readline(io))
        line2 = split(readline(io))
        # Read all lines and flatten the data. This accepts an old superfish format
        dat = readdlm(io)
    end

    # The length of the second line gives a clue about the problem type
    if length(line2) == 1
        problem = "fish"
        line3 = dat[1:3]  # Final header line
        dat = dat[4:end]  # body data
        n = length(dat)
        n % 4 == 0 || error("$n should be divisible by 4")
        dat = reshape(dat, div(n, 4), 4)
    else
        problem = "poisson"
        n = length(dat)
        n % 2 == 0 || error("$n should be divisible by 2")
        dat = reshape(dat, div(n, 2), 2)
    end

    components = Dict{String,Any}()
    if problem == "fish"
        # zmin(cm), zmax(cm), nz-1
        # freq(MHz)
        # rmin(cm), rmax(cm), nr-1
        # 4 columns of data: Ez(MV/m), Er(MV/m), E(MV/m), Hphi(A/m)

        # FISH problem
        zmin, zmax, nz = (
            parse(Float64, line1[1]) * 1e-2,
            parse(Float64, line1[2]) * 1e-2,
            parse(Int, line1[3]) + 1,
        )
        frequency = parse(Float64, line2[1]) * 1e6  # MHz -> Hz
        _rmin, rmax, nr = (
            parse(Float64, line3[1]) * 1e-2,
            parse(Float64, line3[2]) * 1e-2,
            parse(Int, line3[3]) + 1,
        )

        # Read and reshape
        dat = reshape(dat, nr, 1, nz, 4)

        components["electricField/z"] = dat[:, :, :, 1] .* 1e6  # MV/m -> V/m
        components["electricField/r"] = dat[:, :, :, 2] .* 1e6  # MV/m -> V/m
        components["magneticField/theta"] = dat[:, :, :, 4] .* -im .* mu_0  # A/m -> T

    else
        # rmin(cm), rmax(cm), nx-1    # r in cylindrical geometry
        # zmin(cm), zmax(cm), ny-1    # z in cylindrical geometry

        # POISSON problem
        _rmin, rmax, nr = (
            parse(Float64, line1[1]) * 1e-2,
            parse(Float64, line1[2]) * 1e-2,
            parse(Int, line1[3]) + 1,
        )
        zmin, zmax, nz = (
            parse(Float64, line2[1]) * 1e-2,
            parse(Float64, line2[2]) * 1e-2,
            parse(Int, line2[3]) + 1,
        )
        frequency = 0

        # The structure here is different
        dat = reshape(dat, nz, 1, nr, 2)

        # type must be specified
        if type == "electric"
            components["electricField/r"] = permutedims(dat[:, :, :, 1], [3, 2, 1])  # V/m
            components["electricField/z"] = permutedims(dat[:, :, :, 2], [3, 2, 1])  # V/m
        elseif type == "magnetic"
            components["magneticField/r"] = permutedims(dat[:, :, :, 1], [3, 2, 1]) .* 1e-4  # G -> T
            components["magneticField/z"] = permutedims(dat[:, :, :, 2], [3, 2, 1]) .* 1e-4  # G -> T
        else
            error("Poisson problems must specify type as 'electric' or 'magnetic'")
        end
    end

    dz = (zmax - zmin) / (nz - 1)
    dr = rmax / (nr - 1)

    # Attributes
    attrs = Dict{String,Any}(
        "eleAnchorPt" => "beginning",
        "gridGeometry" => "cylindrical",
        "axisLabels" => ("r", "theta", "z"),
        "gridLowerBound" => (0, 1, 0),
        "gridSize" => (nr, 1, nz),
        "gridSpacing" => (dr, 0, dz),
        "gridOriginOffset" => (0, 0, zmin),
        "fundamentalFrequency" => frequency,
        "RFphase" => 0
    )

    if frequency == 0
        attrs["harmonic"] = 0
    else
        attrs["harmonic"] = 1
    end

    return Dict("attrs" => attrs, "components" => components)
end 