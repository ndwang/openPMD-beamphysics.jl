# Genesis 1.3
# -------------
# Version 2 routines

"""
    genesis2_beam_data1(pg)

Calculate statistics of a single particlegroup,
for use in a Genesis 1.3 v2 beam file

Returns a dict of:
    zpos or tpos      : z or t of the mean in m or s
    curpeak           : current in A
    gamma0            : average relativistic gamma (dimensionless)
    emitx, emity      : Normalized emittances in m*rad
    rxbeam, rybeam    : sigma_x, sigma_y in m
    xbeam, ybeam      : <x>, <y> in m
    pxbeam, pybeam    : beta_x*gamma, beta_y*gamma (dimensionless)
    alpha_x, alpha_y  : Twiss alpha  (dimensionless)
"""
function genesis2_beam_data1(pg)
    d = Dict{String,Any}()

    # Handle z or t
    if length(unique(pg.z)) == 1
        # should be different t
        d["tpos"] = pg["mean_t"]
    elseif length(unique(pg.t)) == 1
        # should be different z
        d["zpos"] = pg["mean_z"]
    else
        error("$(pg) has mixed t and z coordinates.")
    end

    d["curpeak"] = pg["average_current"]
    d["gamma0"] = pg["mean_gamma"]
    d["delgam"] = pg["sigma_gamma"]
    d["emitx"] = pg["norm_emit_x"]
    d["emity"] = pg["norm_emit_y"]
    d["rxbeam"] = pg["sigma_x"]
    d["rybeam"] = pg["sigma_y"]
    d["xbeam"] = pg["mean_x"]
    d["ybeam"] = pg["mean_y"]

    d["pxbeam"] = pg["mean_px"] / pg.mass  # beta_x*gamma
    d["pybeam"] = pg["mean_py"] / pg.mass  # beta_y*gamma

    # Twiss, for alpha
    twiss_x = twiss_calc(pg.cov("x", "xp"))
    twiss_y = twiss_calc(pg.cov("y", "yp"))
    d["alphax"] = twiss_x["alpha"]
    d["alphay"] = twiss_y["alpha"]

    return d
end

"""
    genesis2_beam_data(pg, n_slice=nothing)

Slices a particlegroup into n_slice and forms the beam columns.

n_slice is the number of slices. If not given, the beam will be divided
so there are 100 particles in each slice.

Returns a dict of beam_columns, for use with write_genesis2_beam_file

See: genesis2_beam_data1
"""
function genesis2_beam_data(pg, n_slice=nothing)
    # Handle z or t
    if length(unique(pg.z)) == 1
        # should be different t
        slice_key = "t"
    elseif length(unique(pg.t)) == 1
        # should be different z
        slice_key = "z"
    else
        error("$(pg) has mixed t and z coordinates.")
    end

    # Automatic slicing
    if isnothing(n_slice)
        n_slice = length(pg) ÷ 100
    end

    # Slice (split)
    pglist = split(pg, n_slice, key=slice_key)

    d = Dict{String,Vector{Float64}}()
    # Loop over the slices
    for pg in pglist
        data1 = genesis2_beam_data1(pg)
        for (k, v) in data1
            if !haskey(d, k)
                d[k] = Float64[]
            end
            push!(d[k], v)
        end
    end

    return d
end

"""
    write_genesis2_beam_file(fname, beam_columns; verbose=false)

Writes a Genesis 1.3 v2 beam file, using a dict beam_columns

The header will be written as:
? VERSION=1.0
? SIZE=<length of the columns>
? COLUMNS <list of columns
<data>

This is a copy of the lume-genesis routine:

genesis.writers.write_beam_file
"""
function write_genesis2_beam_file(fname::String, beam_columns::Dict; verbose::Bool=false)
    # Get size
    names = collect(keys(beam_columns))
    size = length(beam_columns[names[1]])
    header = """? VERSION=1.0
? SIZE=$size
? COLUMNS $(join(uppercase.(names), " "))"""

    dat = hcat([beam_columns[name] for name in names]...)

    open(fname, "w") do io
        println(io, header)
        writedlm(io, dat, " ")
    end

    verbose && println("Beam written: $fname")

    return header
end

"""
    genesis2_dpa_to_data(dpa, xlamds, current, zsep=1, species="electron")

Converts Genesis 1.3 v2 dpa data to ParticleGroup data.

The Genesis 1.3 v2 dpa phase space coordinates and units are:
    gamma [1]
    phase [rad]
    x [m]
    y [m]
    px/mc [1]
    py/mc [1]
The definition of the phase is different between the .par and .dpa files.
    .par file: phase = psi  = kw*z + field_phase
    .dpa file: phase = kw*z

Parameters
----------
dpa: Array
    Parsed .dpa file as an array with shape (n_slice, 6, n_particles_per_slice)

xlamds: Float64
    wavelength (m)

zsep: Int
    slice separation interval

current: Array
    current array of length n_slice (A)

species: String, required to be 'electron'

Returns
-------
data: Dict with keys: x, px, y, py, z, pz, weight, species, status
    in the units of the openPMD-beamphysics Julia package:
    m, eV/c, m, eV/c, m, eV/c, C

    These are returned in z-coordinates, with z=0.
"""
function genesis2_dpa_to_data(dpa::Array, xlamds::Float64, current::Array, zsep::Int=1, species::String="electron")
    species == "electron" || error("Only electron species is supported")
    mc2 = mec2

    dz = xlamds * zsep

    nslice, dims, n1 = size(dpa)  # n1 particles in a single slice
    dims == 6 || error("Expected 6 dimensions in dpa data")
    n_particle = n1 * nslice

    gamma = vec(dpa[:, 1, :])
    phase = vec(dpa[:, 2, :])
    x = vec(dpa[:, 3, :])
    y = vec(dpa[:, 4, :])
    px = vec(dpa[:, 5, :]) .* mc2
    py = vec(dpa[:, 6, :]) .* mc2
    pz = sqrt.((gamma.^2 .- 1) .* mc2^2 .- px.^2 .- py.^2)

    i0 = 0:(nslice-1)
    i_slice = repeat(i0, inner=n1)

    # Spread particles out over zsep interval
    z = dz .* (i_slice .+ mod.(phase ./ (2π * zsep), zsep))
    z = vec(z)

    # z-coordinates
    t = -z ./ c_light
    z = fill(0.0, n_particle)

    weight = repeat(current, inner=n1) .* dz ./ c_light ./ n1

    return Dict{String,Any}(
        "t" => t,
        "x" => x,
        "px" => px,
        "y" => y,
        "py" => py,
        "z" => z,
        "pz" => pz,
        "species" => species,
        "weight" => weight,
        "status" => fill(1, n_particle)
    )
end

# -------------
# Version 4 routines

"""
    genesis4_beam_data(pg, n_slice=nothing)

Slices a particlegroup into n_slice and forms the sliced beam data.

n_slice is the number of slices. If not given, the beam will be divided
so there are 100 particles in each slice.

Returns a dict of beam_columns, for use with write_genesis2_beam_file

This uses the same routines as genesis2, with some relabeling
See: genesis2_beam_data1
"""
function genesis4_beam_data(pg, n_slice=nothing)
    # Use the same routine as genesis2
    return genesis2_beam_data(pg, n_slice)
end

"""
    write_genesis4_beam(particle_group, h5_fname, n_slice=nothing; verbose=false, return_input_str=false)

Writes sliced beam data to an HDF5 file.
"""
function write_genesis4_beam(particle_group, h5_fname::String, n_slice=nothing; verbose::Bool=false, return_input_str::Bool=false)
    beam_data, units = genesis4_beam_data(particle_group, n_slice=n_slice)

    h5open(h5_fname, "w") do h5
        for (k, v) in beam_data
            h5[k] = v
            write_unit_h5(h5[k], units[k])
        end
    end

    verbose && println("Genesis4 beam file written: $h5_fname")

    if return_input_str
        data_keys = collect(keys(beam_data))
        lines = genesis4_profile_file_input_str(data_keys, h5_fname)
        lines *= genesis4_beam_input_str(data_keys)
        return lines
    end
end

"""
    _profile_file_lines(label, h5filename, xdata_key, ydata_key; isTime=false, reverse=false)

Helper function for generating profile file lines.
"""
function _profile_file_lines(label::String, h5filename::String, xdata_key::String, ydata_key::String; isTime::Bool=false, reverse::Bool=false)
    lines = """&profile_file
  label = $label
  xdata = $h5filename/$xdata_key
  ydata = $h5filename/$ydata_key"""
    
    if isTime
        lines *= "\n  isTime = T"
    end
    if reverse
        lines *= "\n  reverse = T"
    end
    lines *= "\n&end\n"
    return lines
end

"""
    genesis4_profile_file_input_str(data_keys, h5filename)

Returns an input str suitable for the main Genesis4 input file for profile data.
"""
function genesis4_profile_file_input_str(data_keys::Vector{String}, h5filename::String)
    # Genesis4 does not understand paths
    h5filename = basename(h5filename)

    if "s" in data_keys
        xdata_key = "s"
        isTime = false
        reverse = false
    elseif "t" in data_keys
        xdata_key = "t"
        isTime = true
        reverse = true
    else
        error("no s or t found")
    end

    lines = ""
    for ydata_key in data_keys
        if ydata_key == xdata_key
            continue
        end
        lines *= _profile_file_lines(ydata_key, h5filename, xdata_key, ydata_key, isTime=isTime, reverse=reverse)
    end

    return lines
end

"""
    genesis4_beam_input_str(data_keys)

Returns an input str suitable for the main Genesis4 input file for profile data.
"""
function genesis4_beam_input_str(data_keys::Vector{String})
    lines = ["&beam"]
    for k in data_keys
        if k in ("s", "t")
            continue
        end
        push!(lines, "  $k = @$k")
    end
    push!(lines, "&end")
    return join(lines, "\n")
end

"""
    write_genesis4_distribution(particle_group, h5file; verbose=false)

Corresponds to the `import distribution` section in the Genesis4 manual.

Writes datasets to an h5 file:

Datasets:
    x: horizontal coordinate in meters
    y: vertical coordinate in meters
    xp: px/pz (dimensionless trace space horizontal momentum)
    yp: py/pz (dimensionless trace space vertical momentum)
    t: time in seconds
    p: relativistic gamma*beta (total momentum divided by mc)

If particles are at different z, they will be drifted to the same z,
because the output should have different times.

If any of the weights are different, the bunch will be resampled
to have equal weights.
"""
function write_genesis4_distribution(particle_group, h5file::String; verbose::Bool=false)
    h5open(h5file, "w") do h5
        if length(unique(particle_group.z)) > 1
            verbose && println("Drifting particles to the same z")
            # Work on a copy, because we will drift
            P = copy(particle_group)
            # Drift to z
            drift_to_z!(P)
        else
            P = particle_group
        end

        if length(unique(P.weight)) > 1
            n = length(P)
            verbose && println("Resampling $n weighted particles")
            P = resample(P, n, equal_weights=true)
        end

        for k in ["x", "xp", "y", "yp", "t"]
            h5[k] = P[k]
        end

        # p is really beta*gamma
        h5["p"] = P["p"] / P.mass
    end

    verbose && println("Datasets x, xp, y, yp, t, p written to: $h5file")
end

"""
    genesis4_par_to_data(h5, species="electron", smear=true)

Converts Genesis 4 data from an h5 handle or file to data for openPMD-beamphysics.

Genesis4 datasets in the HDF5 file are named:
    x: x position in meters
    px: gamma * beta_x
    y: y position in meters
    py: gamma * beta_y
    theta: angle within a slice in radians
    gamma: relativistic gamma
    current: Current in a single slice (scalar) in Amps

Parameters
----------
h5: HDF5 handle or String
    HDF5 file handle or path to file

species: String, default="electron"
    Particle species

smear: Bool, default=true
    Genesis4 often samples the beam by skipping slices in a step called 'sample'.
    This will smear out the theta coordinate over these slices, preserving the modulus.

Returns
-------
data: Dict
    Data dictionary for ParticleGroup
"""
function genesis4_par_to_data(h5, species::String="electron", smear::Bool=true)
    # Allow for opening a file
    if isa(h5, String)
        isfile(h5) || error("File does not exist: $h5")
        h5 = h5open(h5, "r")
    end

    species == "electron" || error("Only electrons supported for Genesis4")

    # Scalar arrays
    scalars = [
        "beamletsize",
        "one4one",
        "refposition",
        "slicecount",
        "slicelength",
        "slicespacing",
        "wavelength"
    ]

    # Read scalar data
    scalar_data = Dict{String,Any}()
    for s in scalars
        if haskey(h5, s)
            scalar_data[s] = read(h5[s])
        end
    end

    # Read particle data
    x = read(h5["x"])
    px = read(h5["px"])
    y = read(h5["y"])
    py = read(h5["py"])
    theta = read(h5["theta"])
    gamma = read(h5["gamma"])
    current = read(h5["current"])

    # Calculate momentum components
    p = gamma .* mec2
    pz = sqrt.(p.^2 .- (px.^2 .+ py.^2) .* mec2^2)
    px = px .* mec2
    py = py .* mec2

    # Handle smearing if requested
    if smear && haskey(scalar_data, "slicespacing")
        dz = scalar_data["slicespacing"]
        theta = mod.(theta .+ dz .* (0:length(theta)-1), 2π)
    end

    # Create data dictionary
    data = Dict{String,Any}(
        "x" => x,
        "px" => px,
        "y" => y,
        "py" => py,
        "pz" => pz,
        "t" => theta ./ (2π * scalar_data["wavelength"]),
        "species" => species,
        "weight" => current ./ length(x),
        "status" => fill(1, length(x))
    )

    return data
end 